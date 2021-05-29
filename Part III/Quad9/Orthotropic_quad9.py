from numpy import array, pi, zeros, ix_, linspace, arange, int32, setdiff1d, unique,linalg, sqrt
from quad9 import quad9_post, quad9_nodal_stress_averaging,quad9_line_load, quad9_orthotropic
from gmsh_post import write_node_data2, write_elements_data
import matplotlib.pylab as plt

LINE_ELEMENT = 8
QUAD_ELEMENT = 10

Empotrado = 1
Natural_Boundary = 2
Placa = 3
Extremos = 4

def quad9_partIII(file_name, properties, properties_load,a):

	fid = open(f'{file_name}',"r")

	#--------------------------------------#
			# NODES - XY #
	while True:
		line = fid.readline() 
		if line.find('$Nodes') >=0:
			break
	Nnodes= int(fid.readline())

	xy = zeros([Nnodes,2])
	for i in range(Nnodes):
		line = fid.readline()
		sl = line.split()
		xy[i,0] = float(sl[1])
		xy[i,1] = float(sl[2])

	fid = open(f'{file_name}',"r")

	#---------------------------------------#
				# ELEMENTS - CONECT #
	while True:
		line = fid.readline()
		if line.find('$Elements') >= 0:
			break
	Nelements = int(fid.readline())

	conec = zeros((Nelements,9), dtype= int32)

	fixed_nodes = []
	natural_nodes = []
	extremos = []
	nextremos = 0
	placa = []
	nplaca = 0
	NQ = 0
	Quads = []

	for i in range(Nelements):
		line = fid.readline()
		sl = line.split()
		element_number = int32(sl[0]) - 1
		element_type = int32(sl[1])
		physical_group = int32(sl[3])
		entity_number = int32(sl[4])

		if element_type == LINE_ELEMENT and physical_group == Natural_Boundary:
			n1 = int32(sl[5]) -1
			n2 = int32(sl[6]) -1
			n3 = int32(sl[7]) -1 
			natural_nodes += [[n1,n2,n3]]

		if element_type == LINE_ELEMENT and physical_group == Empotrado:
			n1 = int32(sl[5]) -1
			n2 = int32(sl[6]) -1
			n3 = int32(sl[7]) -1
			fixed_nodes += [n1,n2,n3]

		if element_type == QUAD_ELEMENT and (physical_group == Placa or physical_group == Extremos):
			n0 = int32(sl[5])-1
			n1 = int32(sl[6])-1
			n2 = int32(sl[7])-1
			n3 = int32(sl[8])-1
			n4 = int32(sl[9])-1
			n5 = int32(sl[10])-1
			n6 = int32(sl[11])-1
			n7 = int32(sl[12])-1
			n8 = int32(sl[13])-1

			NQ+=1
			conec[element_number,:] = [n0, n1, n2, n3, n4, n5, n6, n7, n8]		
			Quads.append(element_number)

			if physical_group == Placa:
				nplaca+=1
				placa.append(element_number)

			if physical_group == Extremos:
				nextremos+=1
				extremos.append(element_number)

	#---------------------------------------#
			# #DIRECT STIFFNESS METHOD #

	NDOFs_per_node = 2
	NDOFs = NDOFs_per_node*(Nnodes)
	K = zeros((NDOFs, NDOFs))
	f = zeros((NDOFs,1))

	for e in Quads:
		ni = conec[e,0]
		nj = conec[e,1]
		nk = conec[e,2]
		nl = conec[e,3]
		nm = conec[e,4]
		nn = conec[e,5]
		no = conec[e,6]
		np = conec[e,7]
		nq = conec[e,8]

		xy_e = xy[[ni, nj, nk, nl, nm, nn, no, np, nq],:]

		if e in placa:
			ke, fe = quad9_orthotropic(xy_e, properties[0])

		if e in extremos:
			ke,fe = quad9_orthotropic(xy_e, properties[1])


		d = [2*ni, 2*ni+1 , 2*nj, 2*nj+1, 2*nk, 2*nk+1, 2*nl, 2*nl+1, 2*nm, 2*nm+1, 2*no, 2*no+1, 2*np, 2*np+1, 2*nq, 2*nq+1]

		for i in range(16): 
			p = d[i]
			for j in range(16): 
				q = d[j]
				K[p,q] += ke[i,j]
			f[p] += fe[i]

	fixed_nodes = unique(fixed_nodes)
	natural_nodes = (natural_nodes)
	constrained_DOFs = []
	for n in fixed_nodes:
		constrained_DOFs += [2*n, 2*n+1]
	free_DOFs = arange(NDOFs)
	free_DOFs = setdiff1d(free_DOFs, constrained_DOFs)

	KFF = K[ix_(free_DOFs,free_DOFs)]
	KFC = K[ix_(free_DOFs,constrained_DOFs)]
	KCF = K[ix_(constrained_DOFs,free_DOFs)]
	KCC = K[ix_(constrained_DOFs,constrained_DOFs)]

	# ----------------
	# con Fuerza 1 kN
	for nn in natural_nodes:
		ni = nn[0]
		nj = nn[1]
		nk = nn[2]
		xy_e = xy[[ni,nj,nk],:]

		fe = quad9_line_load(xy_e,properties_load)
		d = [2*ni, 2*ni+1, 2*nj, 2*nj+1, 2*nk, 2*nk+1]

		for i in range(6):
			p = d[i]
			f[p] += fe[i]
	
	# ----------------

	ff = f[free_DOFs]
	fc = f[constrained_DOFs]

	ff = f[free_DOFs]
	fc = f[constrained_DOFs]

	u = zeros((NDOFs,1))
	R = zeros((len(constrained_DOFs),1))

	from numpy.linalg import solve

	u[free_DOFs] = solve(KFF, ff)
	uf = u[free_DOFs]
	uc = u[constrained_DOFs]

	uv = u.reshape([-1,2])

	#Desplazamientos
	nodes = arange(1,Nnodes+1)
	write_node_data2(f"desplazamientos_orthotropic_{a}_{file_name}", nodes, uv[:,0], uv[:,1], "Desplazamientos")

	#Calculo de Tensiones
	average_stress_x = []
	average_stress_y = []
	average_stress_xy = []

	for e in Quads:
		sigma_x = zeros((9,1))
		sigma_y = zeros((9,1))
		sigma_xy = zeros((9,1))

		ni = conec[e,0]
		nj = conec[e,1]
		nk = conec[e,2]
		nl = conec[e,3]
		nm = conec[e,4]
		nn = conec[e,5]
		no = conec[e,6]
		np = conec[e,7]
		nq = conec[e,8]

		xy_e = xy[[ni, nj, nk, nl, nm, nn, no, np, nq, ni],:]
		uv_e = uv[[ni, nj, nk, nl, nm, nn, no, np, nq], :]

		u_e = uv_e.reshape((18,1))

		gauss_rule = [
		(-sqrt(3/5), -sqrt(3/5), 5/9, 5/9),                   #1  #0
		(sqrt(3/5), -sqrt(3/5), 5/9, 5/9),                    #3  #1
		(sqrt(3/5), sqrt(3/5), 5/9, 5/9),                     #9  #2
		(-sqrt(3/5), sqrt(3/5), 5/9, 5/9),                    #7  #3			
		(0.0, -sqrt(3/5), sqrt(40/81), sqrt(40/81)),          #2  #4		
		(sqrt(3/5), 0.0, sqrt(40/81), sqrt(40/81)),           #6  #5
		(0, sqrt(3/5), sqrt(40/81), sqrt(40/81)),             #8  #6			
		(-sqrt(3/5), 0.0,sqrt(40/81), sqrt(40/81),),          #4  #7
		(0.0, 0.0, sqrt(64/81), sqrt(64/81),),                #5  #8	
		]

		epsilon, sigma = quad9_post(xy_e, u_e, properties[0])

		for xi, eta, wi, wj in gauss_rule:
			if xi == -sqrt(3/5) and eta == -sqrt(3/5):
				epsilon,(sigma_x[0],sigma_y[0],sigma_xy[0]) = quad9_post(xy_e,u_e,properties[0])
			
			if xi == sqrt(3/5) and eta == -sqrt(3/5):
				epsilon,(sigma_x[1],sigma_y[1],sigma_xy[1]) = quad9_post(xy_e,u_e,properties[0])
			
			if xi == sqrt(3/5) and eta == sqrt(3/5):
				epsilon,(sigma_x[2],sigma_y[2],sigma_xy[2]) = quad9_post(xy_e,u_e,properties[0])
			
			if xi == -sqrt(3/5) and eta == sqrt(3/5):
				epsilon,(sigma_x[3],sigma_y[3],sigma_xy[3]) = quad9_post(xy_e,u_e,properties[0])
			
			if xi == 0 and eta == -sqrt(3/5):
				epsilon,(sigma_x[4],sigma_y[4],sigma_xy[4]) = quad9_post(xy_e,u_e,properties[0])
			
			if xi == sqrt(3/5) and eta == 0:
				epsilon,(sigma_x[5],sigma_y[5],sigma_xy[5]) = quad9_post(xy_e,u_e,properties[0])
			
			if xi == 0 and eta == sqrt(3/5):
				epsilon,(sigma_x[6],sigma_y[6],sigma_xy[6]) = quad9_post(xy_e,u_e,properties[0])
			
			if xi == -sqrt(3/5) and eta == 0:
				epsilon,(sigma_x[7],sigma_y[7],sigma_xy[7]) = quad9_post(xy_e,u_e,properties[0])
			
			if xi == 0 and eta == 0:
				epsilon,(sigma_x[8],sigma_y[8],sigma_xy[8]) = quad9_post(xy_e,u_e,properties[0])

		sx = quad9_nodal_stress_averaging(xy_e, properties[0], sigma_x)
		sy = quad9_nodal_stress_averaging(xy_e, properties[0], sigma_y)
		sxy = quad9_nodal_stress_averaging(xy_e, properties[0], sigma_xy)
		average_stress_x.append(max(sx)[0])
		average_stress_y.append(max(sy)[0])
		average_stress_xy.append(max(sxy)[0])


	elementos = array(Quads)+1
	write_elements_data(f"sigma_x_average_{a}_{file_name}", elementos, average_stress_x, "sigma x average")
	write_elements_data(f"sigma_y_average_{a}_{file_name}", elementos, average_stress_y, "sigma y average")

	#write_elements_data(f"sigma_y_{file_name}", elementos, sigma_y, "sigma y")
	#write_elements_data(f"sigma_x_{file_name}", elementos, sigma_x, "sigma x")

	final = print(f"Listo! revisar archivos de {file_name}")
	return final


# -------------------------------------------------------- #

properties0 = {}
properties0["E"] = 35e9
properties0["nu"] = 0.4
properties0['alpha'] = 4
properties0["bx"] = 0
properties0["by"] = 0
properties0["t"] = 4e-3

properties1 = {}
properties1["E"] = 35e9
properties1["nu"] = 0.4
properties1['alpha'] = 4
properties1["bx"] = 0
properties1["by"] = 0
properties1["t"] = 5e-3

properties_load = {}
properties_load["t"] = properties1["t"]
properties_load["tx"] = 1e3/(properties_load["t"]*4)
properties_load["ty"] = 1e3/(properties_load["t"]*4)


properties = [properties0, properties1]

file_name = ["025_quad9.msh"]
for i in file_name:
	print("Empezando...")
	print(f'File: {i}')
	quad9_partIII(i,properties,properties_load,4)

	
# -------------------------------------------------------- #
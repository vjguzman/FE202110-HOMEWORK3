from numpy import array, pi, zeros, ix_, linspace, arange, int32, setdiff1d, unique,linalg
from quad4 import quad4, quad4_post
from gmsh_post import write_node_data2, write_elements_data
import matplotlib.pylab as plt

LINE_ELEMENT = 1
TRI_ELEMENT = 2
QUAD_ELEMENT = 3

Empotrado = 1
Natural_Boundary = 2
Placa = 3
Extremos = 4

fid = open(f'placa_2d_mesh.msh',"r")

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

fid = open(f'placa_2d_mesh.msh',"r")
#---------------------------------------#
			# ELEMENTS - CONECT #
while True:
	line = fid.readline()
	if line.find('$Elements') >= 0:
		break
Nelements = int(fid.readline())

conec = zeros((Nelements,4), dtype= int32)

fixed_nodes = []
NQ = 0
Quads = []
for i in range(Nelements):
	line = fid.readline()
	sl = line.split()
	element_number = int32(sl[0]) - 1
	element_type = int32(sl[1])
	physical_group = int32(sl[3])
	entity_number = int32(sl[4])

	if element_type == LINE_ELEMENT and physical_group == Empotrado:
		n1 = int32(sl[5]) -1
		n2 = int32(sl[6]) -1
		fixed_nodes += [n1,n2]

	if element_type == QUAD_ELEMENT and (physical_group == Placa or physical_group == Extremos):
		n0 = int32(sl[5])-1
		n1 = int32(sl[6])-1
		n2 = int32(sl[7])-1
		n3 = int32(sl[8])-1
		NQ+=1
		conec[element_number,:] =[n0, n1,n2,n3]		
		Quads.append(element_number)

print(fixed_nodes)
#---------------------------------------#
g = 9.81
densidad = 1.24e3

properties = {}
properties["E"] = 35000
properties["nu"] = 0.4
properties["bx"] = 0
#properties["by"] = -g*densidad
properties["by"] = 0
properties["t"] = 3e-3

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

	xy_e = xy[[ni, nj, nk, nl],:]

	ke, fe = quad4(xy_e, properties)


	d = [2*ni, 2*ni+1 , 2*nj, 2*nj+1, 2*nk, 2*nk+1, 2*nl, 2*nl+1]

	for i in range(8): 
		p = d[i]
		for j in range(8): 
			q = d[j]
			K[p,q] += ke[i,j]
		f[p] += fe[i]

fixed_nodes = unique(fixed_nodes)
constrained_DOFs = []
for n in fixed_nodes:
	constrained_DOFs += [2*n, 2*n+1]

free_DOFs = arange(NDOFs)
free_DOFs = setdiff1d(free_DOFs, constrained_DOFs)

KFF = K[ix_(free_DOFs,free_DOFs)]
KFC = K[ix_(free_DOFs,constrained_DOFs)]
KCF = K[ix_(constrained_DOFs,free_DOFs)]
KCC = K[ix_(constrained_DOFs,constrained_DOFs)]

for i in range(Nnodes):
	f[i] = 1e3                 #1 KN


ff = f[free_DOFs]
fc = f[constrained_DOFs]

ff = f[free_DOFs]
fc = f[constrained_DOFs]


'''
fe_ff = zeros((len(free_DOFs),1))
for i in range(len(ff)):
	if i%2 == 0:
		fe_ff[i] = (properties['by']*f[i][0])
	else:
		fe_ff[i] = (properties['bx']*f[i][0])

print(len(fe_ff))
'''

u = zeros((NDOFs,1))
R = zeros((len(constrained_DOFs),1))

from numpy.linalg import solve

u[free_DOFs] = solve(KFF, ff)
uf = u[free_DOFs]
uc = u[constrained_DOFs]

print(uf)

uv = u.reshape([-1,2])

'''
plt.matshow(K)
plt.title(f"Colormap of Matrix K")
#plt.show()

plt.figure()
factor = 1e1
plt.plot(xy[:,0] + factor*uv[:,0], xy[:,1] + factor*uv[:,1] , ".")
for e in Quads:
	ni = conec[e,0]
	nj = conec[e,1]
	nk = conec[e,2]
	nl = conec[e,3]
	xy_e = xy[[ni,nj,nk,nl,ni],:] +  factor*uv[[ni,nj,nk,nl,ni],:]

	plt.plot(xy_e[:,0], xy_e[:,1], "k")
plt.axis("equal")
plt.title(f" Deformed Structure (δFE displacements)")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.show()

'''

#Desplazamientos
nodes = arange(1,Nnodes+1)
write_node_data2(f"desplazamientos_placa", nodes, uv[:,0], uv[:,1], "Desplazamientos")

#Calculo de Tensiones
sigma_x = zeros(NQ+1)
sigma_y = zeros(NQ+1)
sigma_xy = zeros(NQ+1)

i = 0
for e in Quads:
	ni = conec[e,0]
	nj = conec[e,1]
	nk = conec[e,2]
	nl = conec[e,3]

	xy_e = xy[[ni, nj, nk,nl],:]
	uv_e = uv[[ni,nj,nk,nl], :]

	u_e = uv_e.reshape((-1))

	epsilon, sigma = quad4_post(xy_e, u_e, properties)
	sigma_x[i] = epsilon[0]
	sigma_y[i] = epsilon[1]
	sigma_xy[i] = epsilon[2]
	i+=1

elementos = array(Quads)+1
write_elements_data(f"sigma_x_placa", elementos, sigma_x, "sigma_x")
write_elements_data(f"sigma_y_placa", elementos, sigma_y, "sigma_y")
write_elements_data(f"sigma_xy_placa", elementos, sigma_xy, "sigma_yx")

print("FIN")

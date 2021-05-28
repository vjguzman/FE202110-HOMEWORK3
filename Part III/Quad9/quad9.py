from numpy import array, sqrt, zeros, ix_
from scipy.linalg import det, inv

def quad9(xy, properties):

	E = properties["E"]
	ν = properties["nu"]
	bx = properties["bx"]
	by = properties["by"]
	t = properties["t"]

	Eσ = E / (1-ν**2) * array(
		[
		[1 , ν , 0       ]       ,
		[ν , 1 , 0       ]       ,
		[0 , 0 , (1-ν)/2 ]
		])

	x0 = xy[0,0]
	x1 = xy[1,0]
	x2 = xy[2,0]
	x3 = xy[3,0]
	x4 = xy[4,0]
	x5 = xy[5,0]
	x6 = xy[6,0]
	x7 = xy[7,0]
	x8 = xy[8,0]

	y0 = xy[0,1]
	y1 = xy[1,1]
	y2 = xy[2,1]
	y3 = xy[3,1]
	y4 = xy[4,1]
	y5 = xy[5,1]
	y6 = xy[6,1]
	y7 = xy[7,1]
	y8 = xy[8,1]

	ke = zeros((18,18))
	fe = zeros((18,1))


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


	for xi, eta, wi, wj in gauss_rule:

		x = x0*(xi-1)*(eta-1)*xi*eta/4. + x1*(xi+1)*(eta-1)*xi*eta/4. + x2*(xi+1)*(eta+1)*xi*eta/4. + x3*(xi-1)*(eta+1)*xi*eta/4. + x4*(1-xi**2)*(eta-1)*eta/2. + x5*(1-eta**2)*(xi+1)*xi/2. + x6*(1-xi**2)*(eta+1)*eta/2. + x7*(1-eta**2)*(xi-1)*xi/2. + x8*(1-eta**2)*(1-xi**2)
		y = y0*(xi-1)*(eta-1)*xi*eta/4. + y1*(xi+1)*(eta-1)*xi*eta/4. + y2*(xi+1)*(eta+1)*xi*eta/4. + y3*(xi-1)*(eta+1)*xi*eta/4. + y4*(1-xi**2)*(eta-1)*eta/2. + y5*(1-eta**2)*(xi+1)*xi/2. + y6*(1-xi**2)*(eta+1)*eta/2. + y7*(1-eta**2)*(xi-1)*xi/2. + y8*(1-eta**2)*(1-xi**2)

		dx_dxi = eta*x0*xi*(eta - 1)/4 + eta*x0*(eta - 1)*(xi - 1)/4 + eta*x1*xi*(eta - 1)/4 + eta*x1*(eta - 1)*(xi + 1)/4 + eta*x2*xi*(eta + 1)/4 + eta*x2*(eta + 1)*(xi + 1)/4 + eta*x3*xi*(eta + 1)/4 + eta*x3*(eta + 1)*(xi - 1)/4 - eta*x4*xi*(eta - 1) - eta*x6*xi*(eta + 1) + x5*xi*(1 - eta**2)/2 + x5*(1 - eta**2)*(xi + 1)/2 + x7*xi*(1 - eta**2)/2 + x7*(1 - eta**2)*(xi - 1)/2 - 2*x8*xi*(1 - eta**2)
		dx_deta = eta*x0*xi*(xi - 1)/4 + eta*x1*xi*(xi + 1)/4 + eta*x2*xi*(xi + 1)/4 + eta*x3*xi*(xi - 1)/4 + eta*x4*(1 - xi**2)/2 - eta*x5*xi*(xi + 1) + eta*x6*(1 - xi**2)/2 - eta*x7*xi*(xi - 1) - 2*eta*x8*(1 - xi**2) + x0*xi*(eta - 1)*(xi - 1)/4 + x1*xi*(eta - 1)*(xi + 1)/4 + x2*xi*(eta + 1)*(xi + 1)/4 + x3*xi*(eta + 1)*(xi - 1)/4 + x4*(1 - xi**2)*(eta - 1)/2 + x6*(1 - xi**2)*(eta + 1)/2
		dy_dxi = eta*xi*y0*(eta - 1)/4 + eta*xi*y1*(eta - 1)/4 + eta*xi*y2*(eta + 1)/4 + eta*xi*y3*(eta + 1)/4 - eta*xi*y4*(eta - 1) - eta*xi*y6*(eta + 1) + eta*y0*(eta - 1)*(xi - 1)/4 + eta*y1*(eta - 1)*(xi + 1)/4 + eta*y2*(eta + 1)*(xi + 1)/4 + eta*y3*(eta + 1)*(xi - 1)/4 + xi*y5*(1 - eta**2)/2 + xi*y7*(1 - eta**2)/2 - 2*xi*y8*(1 - eta**2) + y5*(1 - eta**2)*(xi + 1)/2 + y7*(1 - eta**2)*(xi - 1)/2
		dy_deta = eta*xi*y0*(xi - 1)/4 + eta*xi*y1*(xi + 1)/4 + eta*xi*y2*(xi + 1)/4 + eta*xi*y3*(xi - 1)/4 - eta*xi*y5*(xi + 1) - eta*xi*y7*(xi - 1) + eta*y4*(1 - xi**2)/2 + eta*y6*(1 - xi**2)/2 - 2*eta*y8*(1 - xi**2) + xi*y0*(eta - 1)*(xi - 1)/4 + xi*y1*(eta - 1)*(xi + 1)/4 + xi*y2*(eta + 1)*(xi + 1)/4 + xi*y3*(eta + 1)*(xi - 1)/4 + y4*(1 - xi**2)*(eta - 1)/2 + y6*(1 - xi**2)*(eta + 1)/2

		dN0_dxi= eta*xi*(eta - 1)/4. + eta*(eta - 1)*(xi - 1)/4.
		dN0_deta= eta*xi*(xi - 1)/4. + xi*(eta - 1)*(xi - 1)/4.
		dN1_dxi= eta*xi*(eta - 1)/4. + eta*(eta - 1)*(xi + 1)/4.
		dN1_deta= eta*xi*(xi + 1)/4. + xi*(eta - 1)*(xi + 1)/4.
		dN2_dxi= eta*xi*(eta + 1)/4. + eta*(eta + 1)*(xi + 1)/4.
		dN2_deta= eta*xi*(xi + 1)/4. + xi*(eta + 1)*(xi + 1)/4.
		dN3_dxi= eta*xi*(eta + 1)/4. + eta*(eta + 1)*(xi - 1)/4.
		dN3_deta= eta*xi*(xi - 1)/4. + xi*(eta + 1)*(xi - 1)/4.
		dN4_dxi= -eta*xi*(eta - 1)
		dN4_deta= eta*(1 - xi**2)/2. + (1/2 - xi**2/2)*(eta - 1)
		dN5_dxi= xi*(1 - eta**2)/2. + (1 - eta**2)*(xi/2 + 1/2)
		dN5_deta= -eta*xi*(xi + 1)
		dN6_dxi= -eta*xi*(eta + 1)
		dN6_deta= eta*(1 - xi**2)/2. + (1 - xi**2)*(eta/2 + 1/2)
		dN7_dxi= xi*(1 - eta**2)/2. + (1/2 - eta**2/2)*(xi - 1)
		dN7_deta= -eta*xi*(xi - 1)
		dN8_dxi= -2*xi*(1 - eta**2)
		dN8_deta= -2*eta*(1 - xi**2)

		#print(f"x = {x} y = {y}")
		#print(f"dx_dxi = {dx_dxi}")
		#print(f"dx_deta = {dx_deta}")
		#print(f"dy_dxi = {dy_dxi}")
		#print(f"dy_deta = {dy_deta}")

		J = array([
		[dx_dxi, dx_deta],
		[dy_dxi, dy_deta]
		]).T

		detJ = det(J)

		if detJ <= 0.:
			print(f"FATAL! detJ <= 0...")
			exit(-1)

		Jinv = inv(J)

		# print(f"J = {J}")
		# print(f"detJ = {detJ}")

		dN0_dxy = Jinv@array([ dN0_dxi, dN0_deta ])
		dN1_dxy = Jinv@array([ dN1_dxi, dN1_deta ])
		dN2_dxy = Jinv@array([ dN2_dxi, dN2_deta ])
		dN3_dxy = Jinv@array([ dN3_dxi, dN3_deta ])
		dN4_dxy = Jinv@array([ dN4_dxi, dN4_deta ])
		dN5_dxy = Jinv@array([ dN5_dxi, dN5_deta ])
		dN6_dxy = Jinv@array([ dN6_dxi, dN6_deta ])
		dN7_dxy = Jinv@array([ dN7_dxi, dN7_deta ])
		dN8_dxy = Jinv@array([ dN8_dxi, dN8_deta ])

		# ε = B ue
		B = zeros((3, 18))
		B[0,0] = dN0_dxy[0]
		B[1,1] = dN0_dxy[1]
		B[2,0] = dN0_dxy[1]
		B[2,1] = dN0_dxy[0]

		B[0,2] = dN1_dxy[0]
		B[1,3] = dN1_dxy[1]
		B[2,2] = dN1_dxy[1]
		B[2,3] = dN1_dxy[0]

		B[0,4] = dN2_dxy[0]
		B[1,5] = dN2_dxy[1]
		B[2,4] = dN2_dxy[1]
		B[2,5] = dN2_dxy[0]

		B[0,6] = dN3_dxy[0]
		B[1,7] = dN3_dxy[1]
		B[2,6] = dN3_dxy[1]
		B[2,7] = dN3_dxy[0]

		B[0,8] = dN4_dxy[0]
		B[1,9] = dN4_dxy[1]
		B[2,8] = dN4_dxy[1]
		B[2,9] = dN4_dxy[0]

		B[0,10] = dN5_dxy[0]
		B[1,11] = dN5_dxy[1]
		B[2,10] = dN5_dxy[1]
		B[2,11] = dN5_dxy[0]

		B[0,12] = dN6_dxy[0]
		B[1,13] = dN6_dxy[1]
		B[2,12] = dN6_dxy[1]
		B[2,13] = dN6_dxy[0]

		B[0,14] = dN7_dxy[0]
		B[1,15] = dN7_dxy[1]
		B[2,14] = dN7_dxy[1]
		B[2,15] = dN7_dxy[0]

		B[0,16] = dN8_dxy[0]
		B[1,17] = dN8_dxy[1]
		B[2,16] = dN8_dxy[1]
		B[2,17] = dN8_dxy[0]

		N0 = (xi-1)*(eta-1)*xi*eta/4
		N1 = (xi+1)*(eta-1)*xi*eta/4
		N2 = (xi+1)*(eta+1)*xi*eta/4
		N3 = (xi-1)*(eta+1)*xi*eta/4
		N4 = (1-xi**2)*(eta-1)*eta/2
		N5 = (1-eta**2)*(xi+1)*xi/2
		N6 = (1-xi**2)*(eta+1)*eta/2
		N7 = (1-eta**2)*(xi-1)*xi/2
		N8 = (1-eta**2)*(1-xi**2)

		N = zeros((2, 18))
		N[0,0] = N0
		N[1,0] = 0
		N[0,1] = 0
		N[1,1] = N0

		N[0,2] = N1
		N[1,2] = 0
		N[0,3] = 0
		N[1,3] = N1
				
		N[0,4] = N2
		N[1,4] = 0
		N[0,5] = 0
		N[1,5] = N2
				
		N[0,6] = N3
		N[1,6] = 0
		N[0,7] = 0
		N[1,7] = N3
				
		N[0,8] = N4
		N[1,8] = 0
		N[0,9] = 0
		N[1,9] = N4
				
		N[0,10] = N5
		N[1,10] = 0
		N[0,11] = 0
		N[1,11] = N5
				
		N[0,12] = N6
		N[1,12] = 0
		N[0,13] = 0
		N[1,13] = N6
				
		N[0,14] = N7
		N[1,14] = 0
		N[0,15] = 0
		N[1,15] = N7

		N[0,16] = N8
		N[1,16] = 0
		N[0,17] = 0
		N[1,17] = N8
		
		#print(f"N = {N}")

		x0 = xy[0,0]
		x1 = xy[1,0]
		x2 = xy[2,0]
		x3 = xy[3,0]
		x4 = xy[4,0]
		x5 = xy[5,0]
		x6 = xy[6,0]
		x7 = xy[7,0]
		x8 = xy[8,0]

		y0 = xy[0,1]
		y1 = xy[1,1]
		y2 = xy[2,1]
		y3 = xy[3,1]
		y4 = xy[4,1]
		y5 = xy[5,1]
		y6 = xy[6,1]
		y7 = xy[7,1]
		y8 = xy[8,1]

		ke += t * wi * wj * B.T @ Eσ @ B * detJ

		#fe += t * wi * wj * N.T @ array([bx,by])* detJ


	return ke,fe


def quad9_post(xy, u_e, properties):

	E = properties["E"]
	ν = properties["nu"]
	bx = properties["bx"]
	by = properties["by"]
	t = properties["t"]

	#Podemos pasarle otros valores de xi y eta
	if "xi" in properties:
		xi = properties["xi"]
	else:
		xi = 0.0

	if "eta" in properties:
		eta = properties["eta"]
	else:
		eta = 0.0

	Eσ = E / (1-ν**2) * array(
		[
		[1 , ν , 0       ]       ,
		[ν , 1 , 0       ]       ,
		[0 , 0 , (1-ν)/2 ]
		])

	x0 = xy[0,0]
	x1 = xy[1,0]
	x2 = xy[2,0]
	x3 = xy[3,0]
	x4 = xy[4,0]
	x5 = xy[5,0]
	x6 = xy[6,0]
	x7 = xy[7,0]
	x8 = xy[8,0]

	y0 = xy[0,1]
	y1 = xy[1,1]
	y2 = xy[2,1]
	y3 = xy[3,1]
	y4 = xy[4,1]
	y5 = xy[5,1]
	y6 = xy[6,1]
	y7 = xy[7,1]
	y8 = xy[8,1]


	x = x0*(xi-1)*(eta-1)*xi*eta/4. + x1*(xi+1)*(eta-1)*xi*eta/4. + x2*(xi+1)*(eta+1)*xi*eta/4. + x3*(xi-1)*(eta+1)*xi*eta/4. + x4*(1-xi**2)*(eta-1)*eta/2. + x5*(1-eta**2)*(xi+1)*xi/2. + x6*(1-xi**2)*(eta+1)*eta/2. + x7*(1-eta**2)*(xi-1)*xi/2. + x8*(1-eta**2)*(1-xi**2)
	y = y0*(xi-1)*(eta-1)*xi*eta/4. + y1*(xi+1)*(eta-1)*xi*eta/4. + y2*(xi+1)*(eta+1)*xi*eta/4. + y3*(xi-1)*(eta+1)*xi*eta/4. + y4*(1-xi**2)*(eta-1)*eta/2. + y5*(1-eta**2)*(xi+1)*xi/2. + y6*(1-xi**2)*(eta+1)*eta/2. + y7*(1-eta**2)*(xi-1)*xi/2. + y8*(1-eta**2)*(1-xi**2)

	dx_dxi = eta*x0*xi*(eta - 1)/4 + eta*x0*(eta - 1)*(xi - 1)/4 + eta*x1*xi*(eta - 1)/4 + eta*x1*(eta - 1)*(xi + 1)/4 + eta*x2*xi*(eta + 1)/4 + eta*x2*(eta + 1)*(xi + 1)/4 + eta*x3*xi*(eta + 1)/4 + eta*x3*(eta + 1)*(xi - 1)/4 - eta*x4*xi*(eta - 1) - eta*x6*xi*(eta + 1) + x5*xi*(1 - eta**2)/2 + x5*(1 - eta**2)*(xi + 1)/2 + x7*xi*(1 - eta**2)/2 + x7*(1 - eta**2)*(xi - 1)/2 - 2*x8*xi*(1 - eta**2)
	dx_deta = eta*x0*xi*(xi - 1)/4 + eta*x1*xi*(xi + 1)/4 + eta*x2*xi*(xi + 1)/4 + eta*x3*xi*(xi - 1)/4 + eta*x4*(1 - xi**2)/2 - eta*x5*xi*(xi + 1) + eta*x6*(1 - xi**2)/2 - eta*x7*xi*(xi - 1) - 2*eta*x8*(1 - xi**2) + x0*xi*(eta - 1)*(xi - 1)/4 + x1*xi*(eta - 1)*(xi + 1)/4 + x2*xi*(eta + 1)*(xi + 1)/4 + x3*xi*(eta + 1)*(xi - 1)/4 + x4*(1 - xi**2)*(eta - 1)/2 + x6*(1 - xi**2)*(eta + 1)/2
	dy_dxi = eta*xi*y0*(eta - 1)/4 + eta*xi*y1*(eta - 1)/4 + eta*xi*y2*(eta + 1)/4 + eta*xi*y3*(eta + 1)/4 - eta*xi*y4*(eta - 1) - eta*xi*y6*(eta + 1) + eta*y0*(eta - 1)*(xi - 1)/4 + eta*y1*(eta - 1)*(xi + 1)/4 + eta*y2*(eta + 1)*(xi + 1)/4 + eta*y3*(eta + 1)*(xi - 1)/4 + xi*y5*(1 - eta**2)/2 + xi*y7*(1 - eta**2)/2 - 2*xi*y8*(1 - eta**2) + y5*(1 - eta**2)*(xi + 1)/2 + y7*(1 - eta**2)*(xi - 1)/2
	dy_deta = eta*xi*y0*(xi - 1)/4 + eta*xi*y1*(xi + 1)/4 + eta*xi*y2*(xi + 1)/4 + eta*xi*y3*(xi - 1)/4 - eta*xi*y5*(xi + 1) - eta*xi*y7*(xi - 1) + eta*y4*(1 - xi**2)/2 + eta*y6*(1 - xi**2)/2 - 2*eta*y8*(1 - xi**2) + xi*y0*(eta - 1)*(xi - 1)/4 + xi*y1*(eta - 1)*(xi + 1)/4 + xi*y2*(eta + 1)*(xi + 1)/4 + xi*y3*(eta + 1)*(xi - 1)/4 + y4*(1 - xi**2)*(eta - 1)/2 + y6*(1 - xi**2)*(eta + 1)/2

	dN0_dxi= eta*xi*(eta - 1)/4. + eta*(eta - 1)*(xi - 1)/4.
	dN0_deta= eta*xi*(xi - 1)/4. + xi*(eta - 1)*(xi - 1)/4.
	dN1_dxi= eta*xi*(eta - 1)/4. + eta*(eta - 1)*(xi + 1)/4.
	dN1_deta= eta*xi*(xi + 1)/4. + xi*(eta - 1)*(xi + 1)/4.
	dN2_dxi= eta*xi*(eta + 1)/4. + eta*(eta + 1)*(xi + 1)/4.
	dN2_deta= eta*xi*(xi + 1)/4. + xi*(eta + 1)*(xi + 1)/4.
	dN3_dxi= eta*xi*(eta + 1)/4. + eta*(eta + 1)*(xi - 1)/4.
	dN3_deta= eta*xi*(xi - 1)/4. + xi*(eta + 1)*(xi - 1)/4.
	dN4_dxi= -eta*xi*(eta - 1)
	dN4_deta= eta*(1 - xi**2)/2. + (1/2 - xi**2/2)*(eta - 1)
	dN5_dxi= xi*(1 - eta**2)/2. + (1 - eta**2)*(xi/2 + 1/2)
	dN5_deta= -eta*xi*(xi + 1)
	dN6_dxi= -eta*xi*(eta + 1)
	dN6_deta= eta*(1 - xi**2)/2. + (1 - xi**2)*(eta/2 + 1/2)
	dN7_dxi= xi*(1 - eta**2)/2. + (1/2 - eta**2/2)*(xi - 1)
	dN7_deta= -eta*xi*(xi - 1)
	dN8_dxi= -2*xi*(1 - eta**2)
	dN8_deta= -2*eta*(1 - xi**2)

	# print(f"x = {x} y = {y}")

	J = array([
	[dx_dxi, dx_deta],
	[dy_dxi, dy_deta]
	]).T

	detJ = det(J)

	if detJ <= 0.:
		print(f"FATAL! detJ <= 0...")
		exit(-1)

	Jinv = inv(J)

	# print(f"J = {J}")
	# print(f"detJ = {detJ}")

	dN0_dxy = Jinv@array([ dN0_dxi, dN0_deta ])
	dN1_dxy = Jinv@array([ dN1_dxi, dN1_deta ])
	dN2_dxy = Jinv@array([ dN2_dxi, dN2_deta ])
	dN3_dxy = Jinv@array([ dN3_dxi, dN3_deta ])
	dN4_dxy = Jinv@array([ dN4_dxi, dN4_deta ])
	dN5_dxy = Jinv@array([ dN5_dxi, dN5_deta ])
	dN6_dxy = Jinv@array([ dN6_dxi, dN6_deta ])
	dN7_dxy = Jinv@array([ dN7_dxi, dN7_deta ])
	dN8_dxy = Jinv@array([ dN8_dxi, dN8_deta ])

	# ε = B ue
	B = zeros((3, 18))
	B[0,0] = dN0_dxy[0]
	B[1,1] = dN0_dxy[1]
	B[2,0] = dN0_dxy[1]
	B[2,1] = dN0_dxy[0]

	B[0,2] = dN1_dxy[0]
	B[1,3] = dN1_dxy[1]
	B[2,2] = dN1_dxy[1]
	B[2,3] = dN1_dxy[0]

	B[0,4] = dN2_dxy[0]
	B[1,5] = dN2_dxy[1]
	B[2,4] = dN2_dxy[1]
	B[2,5] = dN2_dxy[0]

	B[0,6] = dN3_dxy[0]
	B[1,7] = dN3_dxy[1]
	B[2,6] = dN3_dxy[1]
	B[2,7] = dN3_dxy[0]

	B[0,8] = dN4_dxy[0]
	B[1,9] = dN4_dxy[1]
	B[2,8] = dN4_dxy[1]
	B[2,9] = dN4_dxy[0]

	B[0,10] = dN5_dxy[0]
	B[1,11] = dN5_dxy[1]
	B[2,10] = dN5_dxy[1]
	B[2,11] = dN5_dxy[0]

	B[0,12] = dN6_dxy[0]
	B[1,13] = dN6_dxy[1]
	B[2,12] = dN6_dxy[1]
	B[2,13] = dN6_dxy[0]

	B[0,14] = dN7_dxy[0]
	B[1,15] = dN7_dxy[1]
	B[2,14] = dN7_dxy[1]
	B[2,15] = dN7_dxy[0]

	B[0,16] = dN8_dxy[0]
	B[1,17] = dN8_dxy[1]
	B[2,16] = dN8_dxy[1]
	B[2,17] = dN8_dxy[0]


	ε = B @ u_e
	σ = Eσ @ ε

	return ε, σ


def quad9_nodal_stress_averaging(xy, u_e, properties, W):
	E = properties["E"]
	ν = properties["nu"]
	bx = properties["bx"]
	by = properties["by"]
	t = properties["t"]

	#Podemos pasarle otros valores de xi y eta
	if "xi" in properties:
		xi = properties["xi"]
	else:
		xi = 0.0

	if "eta" in properties:
		eta = properties["eta"]
	else:
		eta = 0.0

	Eσ = E / (1-ν**2) * array(
		[
		[1 , ν , 0       ]       ,
		[ν , 1 , 0       ]       ,
		[0 , 0 , (1-ν)/2 ]
		])

	x0 = xy[0,0]
	x1 = xy[1,0]
	x2 = xy[2,0]

	y0 = xy[0,1]
	y1 = xy[1,1]
	y2 = xy[2,1]



	gauss_rule = [
		(-sqrt(3/5), -sqrt(3/5), 5/9, 5/9),                   #1  #0
		(sqrt(3/5), -sqrt(3/5), 5/9, 5/9),                    #3  #1
		(sqrt(3/5), sqrt(3/5), 5/9, 5/9),                     #9  #2
		(-sqrt(3/5), sqrt(3/5), 5/9, 5/9),                    #7  #3			
		(0.0, -sqrt(3/5), sqrt(40/81), sqrt(40/81)),          #2  #4		
		(sqrt(3/5), 0.0, sqrt(40/81), sqrt(40/81)),           #6  #5
		(0, sqrt(3/5), sqrt(40/81), sqrt(40/81)),             #8  #6			
		(-sqrt(3/5), 0.0,sqrt(40/81), sqrt(40/81),),          #4  #7
		(0.0, 0.0, sqrt(64/81), sqrt(64/81),),                 #5  #8	
	]


	for xi, eta, wi, wj in gauss_rule:
		
		xi = wi
		eta = wj

		N0 = (xi-1)*(eta-1)*xi*eta/4
		N1 = (xi+1)*(eta-1)*xi*eta/4
		N2 = (xi+1)*(eta+1)*xi*eta/4
		N3 = (xi-1)*(eta+1)*xi*eta/4
		N4 = (1-xi**2)*(eta-1)*eta/2
		N5 = (1-eta**2)*(xi+1)*xi/2
		N6 = (1-xi**2)*(eta+1)*eta/2
		N7 = (1-eta**2)*(xi-1)*xi/2
		N8 = (1-eta**2)*(1-xi**2)
		
		N = zeros((4, 4))
		N[0,0] = N2
		N[0,1] = N3
		N[0,2] = N0
		N[0,3] = N1
		N[1,0] = N1
		N[1,1] = N2
		N[1,2] = N3
		N[1,3] = N0
		N[2,0] = N0
		N[2,1] = N1
		N[2,2] = N2
		N[2,3] = N3
		N[3,0] = N3
		N[3,1] = N0
		N[3,2] = N1
		N[3,3] = N2


		w = N.T @ W

	return w


def quad9_line_load(xy, properties_load):
	t = properties_load["t"]
	tx = properties_load["tx"]
	ty = properties_load["ty"]

	x0 = xy[0,0]
	x1 = xy[1,0]
	x2 = xy[2,0]

	y0 = xy[0,1]
	y1 = xy[1,1]
	y2 = xy[2,1]

	xi = -(5/9)/sqrt(3/5)
	eta =  -(5/9)/sqrt(3/5)

	N0 = (xi-1)*(eta-1)*xi*eta/4
	N1 = (xi+1)*(eta-1)*xi*eta/4
	N2 = (xi+1)*(eta+1)*xi*eta/4
	N3 = (xi-1)*(eta+1)*xi*eta/4
	N4 = (1-xi**2)*(eta-1)*eta/2
	N5 = (1-eta**2)*(xi+1)*xi/2
	N6 = (1-xi**2)*(eta+1)*eta/2
	N7 = (1-eta**2)*(xi-1)*xi/2
	N8 = (1-eta**2)*(1-xi**2)

	N = zeros((2, 4))
	N[0,0] = 1
	N[0,1] = 0
	N[0,2] = 1
	N[0,3] = 0

	N[1,0] = 0
	N[1,1] = 1
	N[1,2] = 0
	N[1,3] = 1

	# dN1 en dchi es 1 y para dN0 es 1		

	x = N0*x0 + N1*x1 + N2*x2 
	y = N0*y0 + N1*y1 + N2*y2 

	Lx = abs(x0-x1)
	Ly = abs(y0-y1)

	L = Lx+Ly

	Fe = t*(L)* N.T * [tx,ty]

	fe = [Fe[0,0],Fe[0,1],Fe[1,0],Fe[1,1],Fe[2,0],Fe[2,1],Fe[3,0],Fe[3,1]]

	return fe





'''
xy = array([
[0,0],
[2,0],
[2,2],
[0,2],
[1,0],
[2,1],
[1,2],
[0,1],
[1,1], 	])

conec = array([
    [0,4],      
    [4,1],      
    [1,5],      
    [5,2],      
    [2,6],      
    [6,3],      
    [3,7],     
    [7,0]], dtype = int)

properties = {}
properties["E"] = 1.
properties["nu"] = 0.25
properties["bx"] = 0
properties["by"] = 1.
properties["t"] = 1.

ke, fe = quad9(xy, properties)

print(ke)

import matplotlib.pylab as plt
plt.figure()
plt.plot(xy[:,0], xy[:,1], "*")
Nnodes = len(xy)
for i in range(Nnodes):
    plt.text(xy[i,0],xy[i,1]," {} ".format(i),color='black')

Nelements = conec.shape[0]
for i in range(Nelements):
    ni = conec[i,0]
    nj = conec[i,1]

    plt.plot(xy[[ni,nj],0],xy[[ni,nj],1],"--", color='blue')

plt.show()
'''
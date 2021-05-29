import matplotlib.pylab as plt



def max_stress(file_name,msh):
	
	Max_Stress = []
	stress = []

	fid = open(f'{file_name}',"r")
	while True:
		line = fid.readline() 
		if line.find(msh) >=0:
			break
	for i in range(int(msh)):
		line = fid.readline()
		sl = line.split()
		s = float(sl[1])
		stress.append((s))

	Max_Stress.append(max((stress)))

	return Max_Stress
	

average_stress_x = ["sigma_x_average_025_quad4.msh", "sigma_x_average_075_quad4.msh", "sigma_x_average_125_quad4.msh"]
average_stress_y = ["sigma_y_average_025_quad4.msh", "sigma_y_average_075_quad4.msh", "sigma_y_average_125_quad4.msh"]

msh = ["2319","381","177"]

m_averagex = []
m_averagey = []

for i in range(3):
	#normal = max_stress(normal_stress[i],msh[i])
	average_x = max_stress(average_stress_x[i],msh[i])
	#m_normal.append(normal[0])
	m_averagex.append(average_x[0])
	#normal = max_stress(normal_stress[i],msh[i])
	average_y = max_stress(average_stress_y[i],msh[i])
	#m_normal.append(normal[0])
	m_averagey.append(average_y[0])

h = [0.25,0.75,1.25]
plt.figure()
plt.plot(h, m_averagex,linewidth = 1 ,color = "green", label = "Stress X")
plt.plot(h, m_averagey, linewidth = 1 ,color = "red", label = "Stress Y")
plt.xlabel('Mesh size h')
plt.ylabel('Maximum absolute stress components [Pa]')
plt.title('Quad 4')
for i in range(3):
	plt.text(h[i], m_averagex[i],"{}".format(h[i]),color='black',fontsize=8)
	plt.text(h[i], m_averagey[i],"{}".format(h[i]),color='black',fontsize=8)
plt.show()
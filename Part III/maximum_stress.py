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
		stress.append(abs(s))

	Max_Stress.append(max((stress)))

	return Max_Stress
	

normal_stress = ["sigma_x_025_quad9.msh", "sigma_x_075_quad9.msh", "sigma_x_125_quad9.msh"]
average_stress = ["sigma_x_average_025_quad9.msh", "sigma_x_average_075_quad9.msh", "sigma_x_average_125_quad9.msh"]

msh = ["2319","381","177"]

m_normal = []
m_average = []

for i in range(3):
	normal = max_stress(normal_stress[i],msh[i])
	average = max_stress(average_stress[i],msh[i])
	m_normal.append(normal[0])
	m_average.append(average[0])

h = [0.25,0.75,1.25]

plt.figure()
plt.plot(h, m_average, "*",linewidth = 1 ,color = "green", label = "Average Stress")
plt.plot(h, m_average, "--",linewidth = 1 ,color = "green")
plt.xlabel('Mesh size h')
plt.ylabel('Maximum absolute stress components')
plt.title('Quad 9')
for i in range(3):
	plt.text(h[i], m_average[i],"  mesh {}  ".format(i),color='black')
plt.show()
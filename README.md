# FE202110-HOMEWORK3 :smile:
General Objective is to study quadrilateral elements (Quad4 and Quad9), nodal stress-averaging and orthotropic material models. 

![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Imagenes/Placa_Enunciado.png)

## :computer: Geometry of the base-plate
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Imagenes/Geometria_Placa_2D.png)

- Include physical groups for the fixed displacement condition, prescribed traction condition and the two element thicknesses.

## :computer: Mesh  of the base-plate
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Imagenes/Mesh_Placa_2D.png)

- Mesh must consist of first-order quadrilateral elements only .
- Mesh size should be about 0.5mm in the vicinity of the hole and about 2mm elsewhere. 

## :books: Analysis

#### Perform FEM analysis on your mesh.
In this part of the Homework we use a mesh of quadrilateral elements with a vicinity about 0.5 mm in the hole and 2mm elsewhere, because we want it to be more specific with better results in the hole that our structure has. We notice that as we use a 2D geometrical blase-plate, we were able to use the same method as in CST element that we used in Homework 2 but for Quad4 based on the Chapter 17 of Felippa's book.   <br>

So we use the FEM recipe seen in class for each Quad element in the structure.
<br>
<br>
<img src="https://latex.codecogs.com/gif.latex?k^{e}&plus;=w_{i}w_{j}\cdot&space;B^{T}E_{?}B\cdot&space;det(J)" title="k^{e}+=w_{i}w_{j}\cdot B^{T}E_{?}B\cdot det(J)" /><br>
<img src="https://latex.codecogs.com/gif.latex?f^{e}&plus;=w_{i}w_{j}\cdot&space;N^{T}b\cdot&space;det(J)" title="f^{e}+=w_{i}w_{j}\cdot N^{T}b\cdot det(J)" /> <br>
<br>
Where J is the Jacobian is define by
<br>
<br>
<img src="https://latex.codecogs.com/gif.latex?J=&space;\begin{bmatrix}&space;\frac{\partial&space;x}{\partial&space;\xi&space;}&space;&&space;\frac{\partial&space;x}{\partial&space;\eta&space;}\\&space;\frac{\partial&space;y}{\partial&space;\xi&space;}&space;&&space;\frac{\partial&space;y}{\partial&space;\eta&space;}&space;\end{bmatrix}" title="J= \begin{bmatrix} \frac{\partial x}{\partial \xi } & \frac{\partial x}{\partial \eta }\\ \frac{\partial y}{\partial \xi } & \frac{\partial y}{\partial \eta } \end{bmatrix}" /> 
<br>
<br>
where
<br>
<br>
<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;x}{\partial&space;\xi&space;}&space;=&space;\sum_{i=0}^{3}\frac{\partial&space;N^{i}}{\partial&space;\xi&space;}\cdot&space;x^{i}" title="\frac{\partial x}{\partial \xi } = \sum_{i=0}^{3}\frac{\partial N^{i}}{\partial \xi }\cdot x^{i}" />

<br>
To complete the analysis, the quad4.py programming was carried out in classes where gauss quadrature is applied, which is a method of approximation of a definite integral. Where B is the constant deformation-displacement matrix and b is the vector with the internal forces of the body, while wi and wj come from the gauss quadrature obtain by applying the 1D Gauss rules to each independent variable in turn.
<br>

#### Visualize displacements and stress in the center of each element (one per element as in the CST element). 
To be able to visualize the displacements and stresses in the elements, we use the code present in the gmsh_post to write the displacement data obtain when we applied the direct stiffness method and for the compute of the stresses we used the function quad4_post that can be found in the quad4.py 
<br>
We calculate the results for displacement ans stress for each element using the force applied in the structure by the self weight in the y direction (-density per gravity) says results can be found in the carpet "Resultados con Peso Propio".

#### Use the readme to display the deformed shape and stress components for a force F = 1 KN.  
Using the same method we can display the deformed shape and stresses for the structure with a force of 1 kN.
<br>
##### DISPLACEMENTS
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/con%20F%20%3D%201%20kN/Desplazamientos.png)
##### STRESS IN X
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/con%20F%20%3D%201%20kN/Sigma_X.png)
##### STRESS IN XY
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/con%20F%20%3D%201%20kN/Sigma_XY.png)
##### STRESS IN Y
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/con%20F%20%3D%201%20kN/Sigma_Y.png)

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
In this part of the Homework we use a mesh of quadrilateral elements with a vicinity about 0.5 mm in the hole and 2mm elsewhere, because we want it to be more specific with better results in the hole that our structure has. We notice that as we use a 2D geometrical blase-plate, we were able to use the CST element that we used in Homework 2.   <br>

So we use the FEM recipe seen in class. 
<br>
<br>
<img src="https://latex.codecogs.com/gif.latex?k^{e}&plus;=w_{i}w_{j}\cdot&space;B^{T}E_{?}B\cdot&space;det(J)" title="k^{e}+=w_{i}w_{j}\cdot B^{T}E_{?}B\cdot det(J)" /><br>
<img src="https://latex.codecogs.com/gif.latex?f^{e}&plus;=w_{i}w_{j}\cdot&space;N^{T}b\cdot&space;det(J)" title="f^{e}+=w_{i}w_{j}\cdot N^{T}b\cdot det(J)" /> <br>
<br>
Where J is the Jacobian is define by
<br>
<img src="https://latex.codecogs.com/gif.latex?J=&space;\begin{bmatrix}&space;\frac{\partial&space;x}{\partial&space;\xi&space;}&space;&&space;\frac{\partial&space;x}{\partial&space;\eta&space;}\\&space;\frac{\partial&space;y}{\partial&space;\xi&space;}&space;&&space;\frac{\partial&space;y}{\partial&space;\eta&space;}&space;\end{bmatrix}" title="J= \begin{bmatrix} \frac{\partial x}{\partial \xi } & \frac{\partial x}{\partial \eta }\\ \frac{\partial y}{\partial \xi } & \frac{\partial y}{\partial \eta } \end{bmatrix}" />          where          <img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;x}{\partial&space;\xi&space;}&space;=&space;\sum_{i=0}^{3}\frac{\partial&space;N^{i}}{\partial&space;\xi&space;}\cdot&space;x^{i}" title="\frac{\partial x}{\partial \xi } = \sum_{i=0}^{3}\frac{\partial N^{i}}{\partial \xi }\cdot x^{i}" />



#### Visualize displacements and stress in the center of each element (one per element as in the CST element). 
#### Use the readme to display the deformed shape and stress components for a force F = 1 KN.  

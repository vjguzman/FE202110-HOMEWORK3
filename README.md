# FE202110-HOMEWORK3 :smile:
General Objective is to study quadrilateral elements (Quad4 and Quad9), nodal stress-averaging and orthotropic material models. 

![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Geometria/Placa_Enunciado.png)

### Assigments
- [x] Part I. Study geometry of the 2-D base plate performing FEM analysis.
- [x] Part II. Study the convergence rate of nodal-averaged stresses for this problem. 
- [ ] Part III. Implement and study the convergence rate of higher order elements and stress fields.


## :computer: Geometry of the base-plate
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Geometria/Geometria_Placa_2D.png)

- Include physical groups for the fixed displacement condition, prescribed traction condition and the two element thicknesses.

## :computer: Mesh  of the base-plate
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Geometria/Mesh_Placa_2D.png)

- Mesh must consist of first-order quadrilateral elements only .
- Mesh size should be about 0.5mm in the vicinity of the hole and about 2mm elsewhere. 

<br>
<br>

## :books: Analysis - Part I

#### ▪ Perform FEM analysis on your mesh.
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

#### ▪ Visualize displacements and stress in the center of each element (one per element as in the CST element). 
To be able to visualize the displacements and stresses in the elements, we use the code present in the gmsh_post to write the displacement data obtain when we applied the direct stiffness method and for the compute of the stresses we used the function quad4_post that can be found in the quad4.py 
<br>


#### ▪ Use the readme to display the deformed shape and stress components for a force F = 1 KN.  
Using the same method we can display the deformed shape and stresses for the structure with a force of 1 kN in the x-direction.
<br>
##### 🔅 DISPLACEMENTS
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Part%20I/Desplazamientos_Malla1.png)
##### 🔅 STRESS IN X
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Part%20I/Sigma_X_Malla1.png)

<br>
<br>

## :books: Analysis - Part II

#### ▪ Stress Recovery 
The elemental nodal point stresses are located on the corners and possibly midpoints of the element, they can be compute at the nodal point from adjacent elements will not generally be the same, because they are not required to be continuous in displacements. <br>
From Chapter 28 of Felippa's Book we know that are two approaches for the to study the recovery of stress values for 2D plane-stress elements. They are: 
<br>
<br>
1️⃣ *Direct evaluation of the stresses*, where we substitute the natural coordinates of the nodal points as       arguments to the shape function modules. The equations used are presented below. <br>
    **Strains** <img src="https://latex.codecogs.com/gif.latex?e&space;=&space;B\cdot&space;u^{e}" title="e = B\cdot u^{e}" /> <br>
    **Stresses** <img src="https://latex.codecogs.com/gif.latex?\sigma&space;=&space;E\cdot&space;e&space;=&space;E\cdot&space;B\cdot&space;u^{e}" title="\sigma = E\cdot e = E\cdot B\cdot u^{e}" />
    <br>
   Where u_e is the vector of computed element node displacement and B is the strain-displacement matrix assembled        with the x and y derivatives of the element shape functions evaluated at the point where the element is. 
<br>
<br>
2️⃣ *Gauss integration points*, where we use the element stiffness integration rule and then extrapolate to the element node points. This approach gives better stress values for quadrilateral elements, so we use this method.
<br>

#### ▪ Implementation of nodal stress averaging 
To understand the extrapolation procedure, it is convenient to consider that the region enclosed by the Gaussian points are the internal elements. So the Gauss element (e ’) is a quadrilateral with four nodes, where its natural quadrilateral coordinates are denoted ξ’ and η ’. Which are related to ξ and η with the following relations.
<br>
So we create a function called **nodal_stress_averaging** in out *quad4.py* file, where the inputs are xy, u_e, properties and the vector of the stresses calculated using the direct evaluation. <br>
Here we define the Gaussian rule with respect to what is present in the figure below, then for each value of ξ and η  we calculate the Ni to build the vector N, so we can calculate the average nodal stress. Finally, we only write the data from the stress in the x direction using the *gmsh_post.py* file. <br>

 ![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Geometria/Gauss_Points_Equation.png)
<br>
◽  **N vector** 
<br>
<img src="https://latex.codecogs.com/gif.latex?N&space;=&space;\begin{bmatrix}&space;N2&space;&&space;N3&space;&&space;N0&space;&&space;N1\\&space;N1&space;&&space;N2&space;&&space;N3&space;&&space;N0\\&space;N0&space;&&space;N1&space;&&space;N2&space;&&space;N3&space;\\&space;N3&space;&&space;N0&space;&&space;N1&space;&&space;N2&space;\end{bmatrix}" title="N = \begin{bmatrix} N2 & N3 & N0 & N1\\ N1 & N2 & N3 & N0\\ N0 & N1 & N2 & N3 \\ N3 & N0 & N1 & N2 \end{bmatrix}" />
<br>
<br>
◽    **Gaussian rule**  
  ![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Geometria/Gauss.png)
<br>
<br>
<br>

#### ▪ Discuss of the results
Note that the same applied load of 1 kN of part I was used.

##### 🔅 STRESS IN X - SIMPLE MESH
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Part%20II/Imagenes/Sigma_X_Simple.png)

<br>
<br>

##### 🔅 STRESS IN Y - SIMPLE MESH
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Part%20II/Imagenes/Sigma_Y_Simple.png)

<br>
<br>

##### 🔅 STRESS IN X - MEDIUM MESH
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Part%20II/Imagenes/Sigma_X_Media.png)

<br>
<br>

##### 🔅 STRESS IN Y - MEDIUM MESH
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Part%20II/Imagenes/Sigma_Y_Media.png)

<br>
<br>

##### 🔅 STRESS IN X - FINE MESH
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Part%20II/Imagenes/Sigma_X_Fina.png)

<br>
<br>

##### 🔅 STRESS IN Y - FINE MESH
![img](https://github.com/vjguzman/FE202110-HOMEWORK3/blob/main/Part%20II/Imagenes/Sigma_Y_Fina.png)

<br>
<br>

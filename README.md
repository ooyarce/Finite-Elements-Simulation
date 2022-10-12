# Homework-3---Part-3

## A) Using the same geometry as in the rest of homework 3.

### Quad9 implementation:

![image](https://user-images.githubusercontent.com/69252038/120533685-007fca00-c3af-11eb-9d85-e687c067546e.png)

- We can see that the implementation doesn´t affect too much in the distribution of the stresses. But this helps to define in a better way the perimeter of the hole in the beam. 

### Coarse, medium and fine meshes for the Quad4 and Quad9 discretizations:

#### Coarse meshes:

- Quad9:

![image](https://user-images.githubusercontent.com/69252038/120534793-3e312280-c3b0-11eb-9f14-2651726e39f6.png)

- Quad4:

![image](https://user-images.githubusercontent.com/69252038/120535810-6b320500-c3b1-11eb-9d87-73b0b83dd951.png)

#### Medium meshes:

- Quad9:

![image](https://user-images.githubusercontent.com/69252038/120536171-d67bd700-c3b1-11eb-96b3-71ec2e29df54.png)

- Quad4:

![image](https://user-images.githubusercontent.com/69252038/120536447-22c71700-c3b2-11eb-8aa9-2dca749c6d83.png)

#### Fine meshes:

- Quad9:

![image](https://user-images.githubusercontent.com/69252038/120536756-776a9200-c3b2-11eb-8afb-ce097129ac75.png)

- Quad4:

![image](https://user-images.githubusercontent.com/69252038/120536992-c0bae180-c3b2-11eb-927f-d7441d67ac06.png)

### Discussion
As we can see, the implementation of the quad9 gives more reliable information of the stresses that suffer the beam. This is seen in the colors of the different graphics, in meshes more coarse and with quad4 its difficult to identify some zones of the beam that have an homogeneous, but in meshes with quad9 this is easy to see. 


## B) Implement Quad4 and Quad9 elements that accept an orthotropic material model (assume EL=αET , with α={1,2,4} )

We take the fine mesh in both cases to notice more easily the differences.

####  α = 1

- Quad9:

![q9_fine_alpha1](https://user-images.githubusercontent.com/69252038/120539141-50fa2600-c3b5-11eb-91c9-14eda0f6dae6.png)


- Quad4:

![q4_fine_alpha1](https://user-images.githubusercontent.com/69252038/120539248-6d965e00-c3b5-11eb-8a36-14afab308523.png)


####  α = 2

- Quad9:

![q9_fine_alpha2](https://user-images.githubusercontent.com/69252038/120539354-80109780-c3b5-11eb-9574-3d442528b971.png)

- Quad4:

![q4_fine_alpha2](https://user-images.githubusercontent.com/69252038/120539383-87d03c00-c3b5-11eb-9460-f70df89ad50a.png)

####  α = 4

- Quad9:

![q9_fine_alpha4](https://user-images.githubusercontent.com/69252038/120539419-90287700-c3b5-11eb-97f7-97ca9f3eea03.png)

- Quad4:

![q4_fine_alpha4](https://user-images.githubusercontent.com/69252038/120539445-97e81b80-c3b5-11eb-8059-0a2c37937948.png)

### Discussion
We can see that the meshes with quad4 have "more resistance" to the same loads apply. This is evident because the quad9 graphics shows more deformation than de quad4 graphics. But it occurs a similar phenomenom, that when we increase the value of α, the beam suffers less deformation and with the same loads apply.


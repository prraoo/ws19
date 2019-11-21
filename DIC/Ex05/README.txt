Assignment H5
TEAM:Pramod; Duy;

Files:
Ex04_Pramod_Duy

3.a) 
The first difference in implementations of ced and ced_fed is in the main() function. In case of
the ced diffusion is carried every iteration given kmax = stopping_time/ time step. In otherwords, 
the image is updated by the function cediff() kmax times while incase of ced_fed, there is no updates. ced_fed takes
stopping_time(T) and cycles (M) and updates the image within the function cediff(). 

The second difference is in the evolution of diffusion image in the function cediff(). In case of
ced_fed, there are two control loops. First loop controls the number of cycles the diffusion process
to be performed and the second loop controls the number of intermediate update steps based on N (min
number of steps based on stopping time, cycles and tau=0.25). In this implementation, the image 
update is done by the outer loop. In contrast, the ced implementation has no explicit control loops 
and the image is updated every step.

3.b)  

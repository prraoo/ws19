Assignment H5
TEAM:Pramod; Duy;

Files:
Ex04_Pramod_Duy
├── 3b_ced_fed_8_1.pgm
├── 3b_ced_fed_8_20.pgm
├── 4b_out_0.1.pgm
├── 4b_out_2.pgm
├── 4b_out_150.pgm
└── isonondiff_decorr.c

Parameters:
> 3b_ced_fed_8_1.pgm
	stopping time T = 8; FED cycles M = 1;
> 3b_ced_fed_8_20.pgm
	stopping time T = 8; FED cycles M = 20;

> 4b_out_0.1.pgm
	lambda = 0.1; sigma = 1.5; step size = 0.1; stopping time = 0.4; correlation = 0.58
> 4b_out_2.pgm
	lambda = 2; sigma = 1.5; step size = 0.1; stopping time = 18.4; correlation = 0.037
> 4b_out_150.pgm
	lambda = 150; sigma = 1.5; step size = 0.1; stopping time = 6.5; correlation = 0.046
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

3.b) On comparing the results of ced with T= 8, and ced_fed with different number of the FED cycles,
the images look very similar.

4.b) 
l=lambda; cov=modulus of the correlation , n=number of iterations

|l	|cor	|n	|
|0.1|0.58	|4	|
|1	|0.035	|384|
|2	|0.037	|184|
|5	|0.041  |94	|
|50	|0.046	|66	|
|150|0.046	|65	|

The table above summarises the changes in image quality and diffusion time. For very low lambda, 
say 0.1, the noise and image features are enhanced together hence a correlation between noise and
image is high. On the other side, if the lambda is very high, for instance 150, both noise and 
image are smoothened hence edges are lost as well.

Only when the lambda is 1 or 2 the noise filtered while preserving the edges of the image. 



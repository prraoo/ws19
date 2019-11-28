Assignment H6
TEAM:Pramod; Duy;

Files:
Ex06_Pramod_Duy/
├── diff_gt_1M.pgm
├── diff_gt_t100.pgm
├── diffreac.c
├── out_2_1M_0.2_100.pgm
├── out_2_20_0.2_100.pgm
└── README.txt
Parameters:

> out_2_20_0.2_100.pgm
lambda=2; alpha=20; time step=0.2; iterations=100

> out_2_1M_0.2_100.pgm
lambda=0.2; alpha=10^6; time step=0.2; iterations=100

> diff_gt_t100.pgm
The difference image between house.pgm and out_2_20_0.2_100.pgm which is result of diffusion 
reaction process. 

> diff_gt_1M.pgm
The difference image between house.pgm and out_2_1M_0.2_100.pgm which is result of pure diffusion 
process.

4.b) The number of interations was chosen such that it is greater than or equal to alpha/time_step

4.c) The parameter alpha was chosen such that it is sufficently a large number such that the 
diffusion-reaction almost becomes zero. In the experiment, alpha was 10^6.

4.d) From the differences we can see that the diffusion-reaction process removes the noise from the
image as well as the smoothens the edges. While the pure-diffusion process removes the sharp edges
and does not remove the noise to a great extent.



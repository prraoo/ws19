Assignment H9
TEAM:Pramod; Duy;

Three steady state images for 3(b):
1) out_138_0.1_10000_1.pgm
	initial image: const138.pgm; step_size = 0.1, iter = 10000, offset = 1
2) out_080_0.1_10000_1.pgm
	initial image: const080.pgm; step_size = 0.1, iter = 10000, offset = 1
3) out_noise_0.1_10000_1.pgm
	initial image: noise.pgm; step_size = 0.1, iter = 10000, offset = 1


The filtered images for 3(c):
1) out_shadow_100.pgm
	initial image: shadow.pgm; step_size = 0.1, iter = 100, offset = 1
2) out_shadow_6E6.pgm
	initial image: shadow.pgm; step_size = 0.1, iter = 6*10^6, offset = 1
	
b) In the compatible case, we have the final image as the multiplicative constant
of the guidance image. The constant is the ratio of average grey  value of 
final image to guidance image. In the results we observe similar results.

The image const080.pgm is brighter than const138.pgm hence the final image
out_138_0.1_10000_1.pgm is brighter than out_080_0.1_10000_1.pgm. The average
pixel value of noise and the const138.pgm are similar hence the final images are
of similar brightness. Hence it shows, irrespective of the inital image, the average
pixels value and the guidance image drives the linear osmosis process.


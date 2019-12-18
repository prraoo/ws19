Assignment H8
TEAM:Pramod; Duy;

Nine filtered versions of kitten.pgm
1. out_kit_100_10_10.pgm:
	tonal weight=100, spatial weight=10, iter=10
2. out_kit_100_10_30.pgm
	tonal weight=100, spatial weight=10, iter=300
3. out_kit_100_1_1.pgm
	tonal weight=100, spatial weight=1, iter=1
4. out_kit_10_1_1.pgm
	tonal weight=100, spatial weight=1, iter=1
5. out_kit_300_10_2.pgm
	tonal weight=300, spatial weight=10, iter=2
6. out_kit_50_10_3.pgm
	tonal weight=50, spatial weight=10, iter=3
7. out_kit_100_10_1.pgm
	tonal weight=100, spatial weight=10, iter=1
8. out_kit_100_10_3.pgm
	tonal weight=100, spatial weight=10, iter=3
9. out_kit_10_10_5.pgm
	tonal weight=10, spatial weight=10, iter=5

b) 
1. The difference between the bright pixels of the whisker and dark pixels around 
them is large, hence the tonal distance is large resulting in very small tonal
weight hence the bright pixels are preserved.

2. By increasing the value of theta_s, we can obtain results similar to 
gaussian filter.

The difference is that the bilateral filter with respect to gaussian filter
takes into account the difference in value with the neighbors to preserve edges 
while smoothing. As theta_s increases, the bilateral filter gradually approximates 
Gaussian convolution more closely because the range of g widens and flattens,the 
intensity of the image.

3. Range filter would be a averaging filter based on the dissimilarity in the
image pixel values. This filter would be similar to isotropic nonlinear diffusion
filter that preserves edges. This can be obtained by increasing theta_t and 
decreasing theta_s

c) 
Filtered version of dog: out_dog_50_2_3.pgm
	tonal weight=50, spatial weight=2, iter=3

d)
Filtered version of kitten: final_kit_50_2_5.pgm
	tonal weight=50, spatial weight=2, iter=5

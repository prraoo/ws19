Assignment H11
TEAM:Pramod; Duy;

4.a) 
Iterations:
As the number of iterations increase, the edges become blurry. For 10 iterations
the edges are sharp and as the we approach 1000 iterations, the image is a lot
blurry. 
eg: 
- iso_10_100_0.1_10.pgm
alpha=10, lambda = 100, time-step=0.1, iter=10
- iso_10_100_0.1_1000.pgm
alpha=10, lambda = 100, time-step=0.1, iter=1000

Alpha:
with small alpha, the flow gradients are almost invisible ( image appears dark) 
wit the increase in alpha, the optic flow is distinctly visible. Without 
the regularization term, the optic flow is not visible
eg: 
- iso_0.001_100_0.1_10.pgm
alpha=0.001, lambda = 100, time-step=0.1, iter=10
- iso_100_100_0.1_10.pgm
alpha=100, lambda = 100, time-step=0.1, iter=10

Lambda: 
For larger lambda, the overall optic flow image gradients are brighter.
eg:
- iso_10_0.1_0.1_10.pgm
alpha=10, lambda = 0.10, time-step=0.1, iter=10
- iso_10_1000_0.1_10.pgm
alpha=10, lambda = 1000, time-step=0.1, iter=10

4.b)
With just the homogeneous regularization, the optic flow image is blurry 
especially around the edges and points like ears of the pig.
eg:
- homo_100_100_0.1_20.pgm
alpha=100, lambda = 100, time-step=0.1, iter=20
- homo_10_100_0.2_10.pgm
alpha=10, lambda = 100, time-step=0.2, iter=10
 


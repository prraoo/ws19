Lambda refers to contrast parameter in the Perona-Malik Diffusivity equation.
Sigma refers to the width of the gaussian kernel used as a regularisation parameter.

+ If the value of lamda is small it means that not only edges of images but also some noise structures will be enhanced. Otherwise, if 
value of lambda is too big, it will lost some information on edges because forward diffusion will be applied on edge whose gradient values is smaller than lambda.

+ If the sigma is too small, then the noise smoothening effect is reduced but if the sigma is too large, then significant 
amount of image information will be lost due to smoothing by the gaussian kernel.

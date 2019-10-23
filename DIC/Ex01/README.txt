Assignment 1

3a) Yes, the evolution of the mean, the maximum, the minimum, and the variance is in accordance with with our expectation. The means remains constant, while min and max move towards the mean with the progress in iterations with reduced variance. 

3b) If the time step is larger than 0.25, it would lead to negative central weights i.e., the stability condition for the central weight is violated. This leads to noisy white images, the results of such a case can be observed in the image: “unstable_progression.png”

Diffusion filtered images:
    • out_ild_t0.2_i40.pgm: Iteration Number= 40; Time step = 0.2
    • out_ild_t0.2_i100.pgm: Iteration Number= 100; Time step = 0.2

Gaussian Convolution Images:
    • out_gauss_sd4_p5.pgm : Standard Deviation Sigma = 4 ; Precision = 5*sigma
    • out_gauss_sd6.325_p5.pgm: Standard Deviation Sigma = 6.325; Precision = 5*sigma

Difference Images:
    • out_diff_i40_p5.pgm
    • out_diff_i100_p5.pgm

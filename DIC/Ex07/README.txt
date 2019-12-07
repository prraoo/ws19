Assignment H7
TEAM:Pramod; Duy;

Files:
├── error-fbs50.pgm
├── error-fista50.pgm
├── error-kacanov.pgm
├── pruebab1-fbs50.pgm
├── pruebab1-fista50.pgm
├── pruebab1-kacanov.pgm
├── pruebab1-ref.pgm
├── README.txt
└── tv-fista

1) pruebab1-ref.pgm vs pruebab1-fbs50.pgm

average absolute difference:   6.5204
maximum absolute difference:  36.0000

2) pruebab1-ref.pgm vs pruebab1-fista50.pgm

average absolute difference:   3.5198
maximum absolute difference:  14.0000

3)pruebab1-ref.pgm vs pruebab1-kacanov.pgm

average absolute difference:   1.8630
maximum absolute difference:  28.0000

From the program difference it is clear that the difference 
between kacanov method and reference FBS method has the least 
difference meaning they are close match. 
Next comes the FISTA algorithm. Finally the largest difference
is between FBS (fbs50) algorithm. 

This indicates that the kacanov method gives the best results. In 
other words, the decreasing order of error can be as follows
FBS(fbs50), FISTA and kacanov.

And between FBS and FISTA we can say that FISTA algorithm gives 
better results in lesser iterations. 

Since L2_Norm(D)^2 is <= 8, this puts a constraint in the step 
size, and the time step is 1/L2_Norm(D)^2.
There limit for the step size, that is 0.125 beyond that the 
algorithm becomes unstable. Whereas the FBS method can have
step size till 0.25.




 

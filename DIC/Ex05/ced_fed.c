#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <limits.h>
#include "FED/fed.h"


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*           COHERENCE-ENHANCING ANISOTROPIC DIFFUSION FILTERING            */
/*                     with FAST EXPLICIT DIFFUSION                         */
/*                                                                          */
/*                 (Copyright by Joachim Weickert 8/2014)                   */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 features:
 - fast explicit diffusion
 - standard discretization of mixed derivatives
 - presmoothing at noise scale:  convolution-based, Neumann b.c.
 - presmoothing at integration scale: convolution-based, Dirichlet b.c.
*/


/*--------------------------------------------------------------------------*/

void alloc_vector

     (float **vector,   /* vector */
      long  n1)         /* size */

     /* allocates memory for a vector of size n1 */


{
*vector = (float *) malloc (n1 * sizeof(float));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough memory available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (float ***matrix,  /* matrix */
      long  n1,         /* size in direction 1 */
      long  n2)         /* size in direction 2 */

     /* allocates memory for matrix of size n1 * n2 */


{
long i;

*matrix = (float **) malloc (n1 * sizeof(float *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough memory available\n");
   exit(1);
   }
for (i=0; i<n1; i++)
    {
    (*matrix)[i] = (float *) malloc (n2 * sizeof(float));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_matrix: not enough memory available\n");
       exit(1);
       }
    }
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

     (float *vector,    /* vector */
      long  n1)         /* size */

     /* disallocates memory for a vector of size n1 */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (float **matrix,   /* matrix */
      long  n1,         /* size in direction 1 */
      long  n2)         /* size in direction 2 */

     /* disallocates memory for matrix of size n1 * n2 */

{
long i;

for (i=0; i<n1; i++)
    free(matrix[i]);

free(matrix);

return;
}

/*--------------------------------------------------------------------------*/

void read_string

     (char *v)         /* string to be read */

/*
 reads a long value v
*/

{
fgets (v, 80, stdin);
if (v[strlen(v)-1] == '\n')
   v[strlen(v)-1] = 0;
return;
}

/*--------------------------------------------------------------------------*/

void read_long

     (long *v)         /* value to be read */

/*
 reads a long value v
*/

{
char   row[80];    /* string for reading data */

fgets (row, 80, stdin);
if (row[strlen(row)-1] == '\n')
   row[strlen(row)-1] = 0;
sscanf(row, "%ld", &*v);
return;
}

/*--------------------------------------------------------------------------*/

void read_float

     (float *v)         /* value to be read */

/*
 reads a float value v
*/

{
char   row[80];    /* string for reading data */

fgets (row, 80, stdin);
if (row[strlen(row)-1] == '\n')
   row[strlen(row)-1] = 0;
sscanf(row, "%f", &*v);
return;
}

/*--------------------------------------------------------------------------*/

void read_pgm_and_allocate_memory

     (const char  *file_name,    /* name of pgm file */ 
      long        *nx,           /* image size in x direction, output */
      long        *ny,           /* image size in y direction, output */
      float       ***u)          /* image, output */   

/* 
  reads a greyscale image that has been encoded in pgm format P5;
  allocates memory for the image u; 
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
FILE   *inimage;    /* input file */
char   row[80];     /* for reading data */
long   i, j;        /* loop variables */

/* open file */
inimage = fopen (file_name, "rb");
if (NULL == inimage) 
   {
   printf ("could not open file '%s' for reading, aborting.\n", file_name);
   exit (1);
   }

/* read header */
fgets (row, 80, inimage);          /* skip format definition */
fgets (row, 80, inimage);        
while (row[0]=='#')                /* skip comments */
      fgets (row, 80, inimage);
sscanf (row, "%ld %ld", nx, ny);   /* read image size */
fgets (row, 80, inimage);          /* read maximum grey value */

/* allocate memory */
alloc_matrix (u, (*nx)+2, (*ny)+2);

/* read image data row by row */
for (j=1; j<=(*ny); j++) 
 for (i=1; i<=(*nx); i++) 
     (*u)[i][j] = (float) getc(inimage);

/* close file */
fclose(inimage);

return;

} /* read_pgm_and_allocate_memory */

/*--------------------------------------------------------------------------*/

void comment_line

     (char* comment,       /* comment string (output) */
      char* lineformat,    /* format string for comment line */
      ...)                 /* optional arguments */

/* 
  Add a line to the comment string comment. The string line can contain plain
  text and format characters that are compatible with sprintf.
  Example call: print_comment_line(comment,"Text %f %d",float_var,int_var);
  If no line break is supplied at the end of the input string, it is added
  automatically.
*/

{
char     line[80];
va_list  arguments;

/* get list of optional function arguments */
va_start(arguments,lineformat);

/* convert format string and arguments to plain text line string */
vsprintf(line,lineformat,arguments);

/* add line to total commentary string */
strncat(comment,line,80);

/* add line break if input string does not end with one */
if (line[strlen(line)-1] != '\n')
   sprintf(comment,"%s\n",comment);

/* close argument list */
va_end(arguments);

return;

} /* comment_line */

/*--------------------------------------------------------------------------*/

void write_pgm

     (float  **u,          /* image, unchanged */ 
      long   nx,           /* image size in x direction */
      long   ny,           /* image size in y direction */
      char   *file_name,   /* name of pgm file */
      char   *comments)    /* comment string (set 0 for no comments) */

/* 
  writes a greyscale image into a pgm P5 file;
*/

{
FILE           *outimage;  /* output file */
long           i, j;       /* loop variables */
float          aux;        /* auxiliary variable */
unsigned char  byte;       /* for data conversion */

/* open file */
outimage = fopen (file_name, "wb");
if (NULL == outimage) 
   {
   printf("Could not open file '%s' for writing, aborting\n", file_name);
   exit(1);
   }

/* write header */
fprintf (outimage, "P5\n");                  /* format */
if (comments != 0)
   fprintf (outimage, comments);             /* comments */
fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");                 /* maximal value */

/* write image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     aux = u[i][j] + 0.499999;    /* for correct rounding */
     if (aux < 0.0)
        byte = (unsigned char)(0.0);
     else if (aux > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(aux);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }

/* close file */
fclose (outimage);

return;

} /* write_pgm */

/*--------------------------------------------------------------------------*/

void dummies
 
     (float **u,        /* image matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

/* creates dummy boundaries by mirroring */

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    u[i][0]    = u[i][1];
    u[i][ny+1] = u[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    u[0][j]    = u[1][j];
    u[nx+1][j] = u[nx][j];
    }
return;
}  

/*--------------------------------------------------------------------------*/

void gauss_conv 

     (float    sigma,     /* standard deviation of Gaussian */
      long     nx,        /* image dimension in x direction */ 
      long     ny,        /* image dimension in y direction */ 
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      float    precision, /* cutoff at precision * sigma */
      long     bc,        /* type of boundary condition */
                          /* 0=Dirichlet, 1=reflecing, 2=periodic */
      float    **f)       /* input: original image ;  output: smoothed */

/* 
 Gaussian convolution.
*/

{
long    i, j, p;              /* loop variables */
long    length;               /* convolution vector: 0..length */
float   sum;                  /* for summing up */
float   *conv;                /* convolution vector */
float   *help;                /* row or column with dummy boundaries */
      

/* ------------------------ diffusion in x direction -------------------- */

/* calculate length of convolution vector */
length = (long)(precision * sigma / hx) + 1;
if ((bc != 0) && (length > nx))
   {
   printf("gauss_conv: sigma too large \n"); 
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length+1);

/* calculate entries of convolution vector */
for (i=0; i<=length; i++)
    conv[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) 
              * exp (- (i * i * hx * hx) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (i=1; i<=length; i++)
    sum = sum + 2.0 * conv[i];
for (i=0; i<=length; i++)
    conv[i] = conv[i] / sum;

/* allocate storage for a row */
alloc_vector (&help, nx+length+length);

for (j=1; j<=ny; j++)
    {
    /* copy in row vector */
    for (i=1; i<=nx; i++)
        help[i+length-1] = f[i][j];

    /* assign boundary conditions */
    if (bc == 0) /* Dirichlet boundary conditions */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = 0.0;
           help[nx+length-1+p] = 0.0;
           }
    else if (bc == 1) /* reflecting b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[length+p-1];
           help[nx+length-1+p] = help[nx+length-p];
           }
    else if (bc == 2) /* periodic b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[nx+length-p];
           help[nx+length-1+p] = help[length+p-1];
           }

    /* convolution step */
    for (i=length; i<=nx+length-1; i++)
        {
        /* calculate convolution */
        sum = conv[0] * help[i];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[i+p] + help[i-p]);
        /* write back */
        f[i-length+1][j] = sum;
        }
    } /* for j */

/* disallocate storage for a row */
disalloc_vector (help, nx+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length + 1);


/* ------------------------ diffusion in y direction -------------------- */

/* calculate length of convolution vector */
length = (long)(precision * sigma / hy) + 1;
if ((bc != 0) && (length > ny))
   {
   printf("gauss_conv: sigma too large \n"); 
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length + 1);

/* calculate entries of convolution vector */
for (j=0; j<=length; j++)
    conv[j] = 1 / (sigma * sqrt(2.0 * 3.1415927)) 
              * exp (- (j * j * hy * hy) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (j=1; j<=length; j++)
    sum = sum + 2.0 * conv[j];
for (j=0; j<=length; j++)
    conv[j] = conv[j] / sum;

/* allocate storage for a row */
alloc_vector (&help, ny+length+length);

for (i=1; i<=nx; i++)
    {
    /* copy in column vector */
    for (j=1; j<=ny; j++)
        help[j+length-1] = f[i][j];

    /* assign boundary conditions */
    if (bc == 0) /* Dirichlet boundary conditions */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = 0.0;
           help[ny+length-1+p] = 0.0;
           }
    else if (bc == 1) /* reflecting b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[length+p-1];
           help[ny+length-1+p] = help[ny+length-p];
           }
    else if (bc == 2) /* periodic b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[ny+length-p];
           help[ny+length-1+p] = help[length+p-1];
           } 
 
    /* convolution step */
    for (j=length; j<=ny+length-1; j++)
        {
        /* calculate convolution */
        sum = conv[0] * help[j];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[j+p] + help[j-p]);
        /* write back */
        f[i][j-length+1] = sum;
        }
    } /* for i */

/* disallocate storage for a row */
disalloc_vector (help, ny+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length+1);

return;

} /* gauss_conv */

/* ------------------------------------------------------------------------ */

void struct_tensor 

     (float    **v,       /* image !! gets smoothed on exit !! */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      float    sigma,     /* noise scale */
      float    rho,       /* integration scale */
      float    **dxx,     /* element of structure tensor, output */
      float    **dxy,     /* element of structure tensor, output */
      float    **dyy)     /* element of structure tensor, output */

/*
 Calculates the structure tensor.
*/

{
long    i, j;                 /* loop variables */
float   dv_dx, dv_dy;         /* derivatives of v */
float   two_hx, two_hy;       /* time savers */


/* ---- smoothing at noise scale, reflecting b.c. ---- */

if (sigma > 0.0) 
   gauss_conv (sigma, nx, ny, hx, hy, 3.0, 1, v);  


/* ---- building tensor product ---- */

two_hx = 2.0 * hx;
two_hy = 2.0 * hy;
dummies (v, nx, ny);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     dv_dx = (v[i+1][j] - v[i-1][j]) / two_hx;
     dv_dy = (v[i][j+1] - v[i][j-1]) / two_hy;
     dxx[i][j] = dv_dx * dv_dx; 
     dxy[i][j] = dv_dx * dv_dy; 
     dyy[i][j] = dv_dy * dv_dy; 
     }


/* ---- smoothing at integration scale, Dirichlet b.c. ---- */

if (rho > 0.0)
   {
   gauss_conv (rho, nx, ny, hx, hy, 3.0, 0, dxx);
   gauss_conv (rho, nx, ny, hx, hy, 3.0, 0, dxy);
   gauss_conv (rho, nx, ny, hx, hy, 3.0, 0, dyy);
   }

return;

} /* struct_tensor */

/* ------------------------------------------------------------------------ */

void PA_trans 

     (float a11,        /* coeffs of (2*2)-matrix */ 
      float a12,        /* coeffs of (2*2)-matrix */ 
      float a22,        /* coeffs of (2*2)-matrix */ 
      float *c,         /* 1. comp. of 1. eigenvector, output */ 
      float *s,         /* 2. comp. of 1. eigenvector, output */ 
      float *lam1,      /* larger  eigenvalue, output */
      float *lam2)      /* smaller eigenvalue, output */

/*
 Principal axis transformation. 
*/

{
float  help, norm;    /* time savers */ 

/* eigenvalues */
help  = sqrt (pow(a22-a11, 2.0) + 4.0 * a12 * a12);
*lam1 = (a11 + a22 + help) / 2.0; 
*lam2 = (a11 + a22 - help) / 2.0; 

/* eigenvectors */
*c = 2.0 * a12;
*s = a22 - a11 + help;

/* normalize eigenvectors */
norm = sqrt (*c * *c + *s * *s);
if (norm >= 0.0000001) 
   {
   *c = *c / norm;   
   *s = *s / norm;   
   }
else  
   {
   *c = 1.0;
   *s = 0.0;
   }
return;

} /* PA_trans */

/* ----------------------------------------------------------------------- */

void PA_backtrans 

     (float  c,      /* 1. comp. of 1. eigenvector */ 
      float  s,      /* 2. comp. of 1. eigenvector */ 
      float  lam1,   /* 1. eigenvalue */ 
      float  lam2,   /* 2. eigenvalue */ 
      float  *a11,   /* coeff. of (2*2)-matrix, output */ 
      float  *a12,   /* coeff. of (2*2)-matrix, output */ 
      float  *a22)   /* coeff. of (2*2)-matrix, output */ 

/*
 Principal axis backtransformation of a symmetric (2*2)-matrix. 
 A = U * diag(lam1, lam2) * U_transpose with U = (v1 | v2)     
 v1 = (c, s) is first eigenvector
*/

{

*a11 = c * c * lam1 + s * s * lam2;
*a22 = lam1 + lam2 - *a11;             /* trace invariance */
*a12 = c * s * (lam1 - lam2);

return;

} /* PA_backtrans */

/*--------------------------------------------------------------------------*/

void diff_tensor 
     
     (float    C,        /* coherence parameter */
      float    alpha,    /* linear diffusion fraction */
      long     nx,       /* image dimension in x direction */
      long     ny,       /* image dimension in y direction */
      float    **dxx,    /* in: structure tensor el., out: diff. tensor el. */
      float    **dxy,    /* in: structure tensor el., out: diff. tensor el. */ 
      float    **dyy)    /* in: structure tensor el., out: diff. tensor el. */ 

/*
 Calculates the diffusion tensor of CED by means of the structure tensor.
*/

{
long    i, j;          /* loop variables */
float   beta;          /* time saver */
float   c, s;          /* specify first eigenvector */
float   mu1, mu2;      /* eigenvalues of structure tensor */
float   lam1, lam2;    /* eigenvalues of diffusion tensor */

beta = 1.0 - alpha;

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* principal axis transformation */
     PA_trans (dxx[i][j], dxy[i][j], dyy[i][j], &c, &s, &mu1, &mu2);

     /* calculate eigenvalues */
     lam1 = alpha;
     lam2 = alpha + beta * exp (- C / pow (mu1-mu2, 2.0));

     /* principal axis backtransformation */
     PA_backtrans (c, s, lam1, lam2, &dxx[i][j], &dxy[i][j], &dyy[i][j]);  
     }

return;

}  /* diff_tensor */

/*--------------------------------------------------------------------------*/

void cediff 

     (long     nx,        /* image dimension in x direction */ 
      long     ny,        /* image dimension in y direction */ 
      float    T,         /* overall stopping time */
      long     M,         /* number of FED cycles */
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      float    C,         /* coherence parameter */
      float    sigma,     /* noise scale */
      float    rho,       /* integration scale */
      float    alpha,     /* linear diffusion fraction */
      float    **u)       /* input: original image;  output: smoothed */

/* 
 Coherence-enhancing anisotropic diffusion. 
 Fast Explicit Diffusion.
*/

{
long    i, j;                 /* loop variables */
float   rxx, rxy, ryy;        /* time savers */
float   wN, wNE, wE, wSE;     /* weights */
float   wS, wSW, wW, wNW;     /* weights */
float   **f;                  /* work copy of u */
float   **dxx, **dxy, **dyy;  /* entries of structure/diffusion tensor */
float   *tau;                 /* Vector of FED time step sizes */
int     N;                    /* Number of steps */
int     n,m;                  /* Loop counters over steps, cycles */


/* ---- allocate storage for f and diffusion tensor entries ---- */

alloc_matrix (&f,   nx+2, ny+2);
alloc_matrix (&dxx, nx+2, ny+2);
alloc_matrix (&dxy, nx+2, ny+2);
alloc_matrix (&dyy, nx+2, ny+2);


/* Initialise step sizes for process with                        */
/* - overall stopping time T                                     */
/* - number of cycles M                                          */
/* - stability limit for 2-D linear diffusion: 0.25              */
N = fed_tau_by_process_time(T, M, 0.25f, 1, &tau);


/* ---- Perform M outer cycles ---- */

  for(m = 0; m < M; m++)
     {
     /* ---- copy u into f ---- */

     for (i=1; i<=nx; i++)
      for (j=1; j<=ny; j++)
          f[i][j] = u[i][j];


     /* ---- calculate entries of structure tensor (alters f!!!) ---- */

     struct_tensor (f, nx, ny, hx, hy, sigma, rho, dxx, dxy, dyy);


     /* ---- calculate entries of diffusion tensor ---- */

     diff_tensor (C, alpha, nx, ny, dxx, dxy, dyy);
     
     /* assign dummy boundaries */
     dummies (dxx, nx, ny);
     dummies (dxy, nx, ny);
     dummies (dyy, nx, ny);


     /* ---- Each cycle performs N steps with varying step size ---- */
     
     for(n = 0; n < N; n++)
        {
        rxx  = tau[n] / (2.0 * hx * hx);
        ryy  = tau[n] / (2.0 * hy * hy);
        rxy  = tau[n] / (4.0 * hx * hy);


        /* copy u into f and assign dummy boundaries */
        for (i=1; i<=nx; i++)
         for (j=1; j<=ny; j++)
             f[i][j] = u[i][j];
        dummies (f, nx, ny);


        /* ---- calculate explicit nonlinear diffusion of u ---- */
        for (i=1; i<=nx; i++)
         for (j=1; j<=ny; j++)
             {
             /* weights */
             wE  =   rxx * (dxx[i+1][j]   + dxx[i][j]);
             wW  =   rxx * (dxx[i-1][j]   + dxx[i][j]);
             wS  =   ryy * (dyy[i][j+1]   + dyy[i][j]);
             wN  =   ryy * (dyy[i][j-1]   + dyy[i][j]);
             wSE =   rxy * (dxy[i+1][j+1] + dxy[i][j]);
             wNW =   rxy * (dxy[i-1][j-1] + dxy[i][j]);
             wNE = - rxy * (dxy[i+1][j-1] + dxy[i][j]);
             wSW = - rxy * (dxy[i-1][j+1] + dxy[i][j]);

             /* modify weights to prevent flux across boundaries */
             if (i==1)  /* set western weights zero */
                {
                wSW = 0.0;
                wW  = 0.0;
                wNW = 0.0;
                }
             if (i==nx) /* set eastern weights zero */
                {
                wNE = 0.0;
                wE  = 0.0;
                wSE = 0.0;
                }
             if (j==1)  /* set northern weights zero */
                {
                wNW = 0.0;
                wN  = 0.0;
                wNE = 0.0;
                }
             if (j==ny) /* set southern weights zero */
                {
                wSE = 0.0;
                wS  = 0.0; 
                wSW = 0.0;
                } 

             /* evolution */
             u[i][j] = f[i][j]  
                     + wE  * (f[i+1][j]   - f[i][j]) 
                     + wW  * (f[i-1][j]   - f[i][j]) 
                     + wS  * (f[i][j+1]   - f[i][j]) 
                     + wN  * (f[i][j-1]   - f[i][j])
                     + wSE * (f[i+1][j+1] - f[i][j]) 
                     + wNW * (f[i-1][j-1] - f[i][j]) 
                     + wSW * (f[i-1][j+1] - f[i][j]) 
                     + wNE * (f[i+1][j-1] - f[i][j]);
             }
        } /* for n steps */
     } /* for m outer cycles*/


/* ---- disallocate storage for f and diffusion tensor entries ---- */

disalloc_matrix (f,   nx+2, ny+2);
disalloc_matrix (dxx, nx+2, ny+2);
disalloc_matrix (dxy, nx+2, ny+2);
disalloc_matrix (dyy, nx+2, ny+2);

printf ("number of iterations: %5ld \n", N*M);

return;

} /* cediff */

/*--------------------------------------------------------------------------*/

void analyse

     (float   **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in y direction */
      float   *min,        /* minimum, output */
      float   *max,        /* maximum, output */
      float   *mean,       /* mean, output */
      float   *std)        /* standard deviation, output */

/*
 computes minimum, maximum, mean, and standard deviation of an image u
*/

{
long    i, j;       /* loop variables */
double  help1;      /* auxiliary variable */
float   help2;      /* auxiliary variable */

*min  = u[1][1];
*max  = u[1][1];
help1 = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     if (u[i][j] < *min) *min = u[i][j];
     if (u[i][j] > *max) *max = u[i][j];
     help1 = help1 + (double)u[i][j];
     }
*mean = (float)help1 / (nx * ny);

*std = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     help2  = u[i][j] - *mean;
     *std = *std + help2 * help2;
     }
*std = sqrt(*std / (nx * ny));

return;

} /* analyse */

/*--------------------------------------------------------------------------*/

int main ()

{
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
float  **u;                  /* image */
long   nx, ny;               /* image size in x, y direction */ 
float  sigma;                /* noise scale */
float  rho;                  /* integration scale */
float  alpha;                /* linear diffusion fraction */
float  C;                    /* coherence parameter */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */
long   M;                    /* number of FED cycles */
float  T;                    /* stopping time */

printf ("\n");
printf ("COHERENCE-ENHANCING ANISOTROPIC DIFFUSION, WITH FED SCHEME\n\n");
printf ("**********************************************************\n\n");
printf ("    Copyright 2014 by Joachim Weickert                      \n");
printf ("    Dept. of Mathematics and Computer Science               \n");
printf ("    Saarland University, Saarbruecken, Germany            \n\n");
printf ("    All rights reserved. Unauthorized usage,                \n");
printf ("    copying, hiring, and selling prohibited.              \n\n");
printf ("    Send bug reports to                                     \n");
printf ("    weickert@mia.uni-saarland.de                          \n\n");
printf ("**********************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                       ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &u);


/* ---- read parameters ---- */

printf ("overall stopping time T (>0) (float):    ");
read_float (&T);

printf ("number of FED cycles M (integer):        ");
read_long (&M);

printf ("output image (pgm):                      ");
read_string (out);
printf ("\n");

/* initializations */
C = 1;
sigma = 0.5;
alpha = 0.001;
rho = 4;


/* ---- process image ---- */

cediff (nx, ny, T, M, 1.0, 1.0, C, sigma, rho, alpha, u);
    
    
/* ---- analyse filtered image ---- */

/* check minimum, maximum, mean, standard deviation */
analyse (u, nx, ny, &min, &max, &mean, &std);
printf ("minimum:           %8.2f \n", min);
printf ("maximum:           %8.2f \n", max);
printf ("mean:              %8.2f \n", mean);
printf ("standard dev.:     %8.2f \n\n", std);


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# coherence-enhancing diffusion, fed scheme\n");
comment_line (comments, "# initial image: %s\n", in);
comment_line (comments, "# C:             %8.4f\n", C);
comment_line (comments, "# sigma:         %8.4f\n", sigma);
comment_line (comments, "# rho:           %8.4f\n", rho);
comment_line (comments, "# alpha:         %8.4f\n", alpha);
comment_line (comments, "# FED cycles:    %8ld\n", M);
comment_line (comments, "# stopping time: %8.4f\n", T);
comment_line (comments, "# min:           %8.4f\n", min);
comment_line (comments, "# max:           %8.4f\n", max);
comment_line (comments, "# mean:          %8.4f\n", mean);
comment_line (comments, "# standard dev.: %8.4f\n", std);

/* write image */
write_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (u, nx+2, ny+2);

return(0);
}

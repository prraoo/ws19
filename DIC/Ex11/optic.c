#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                               OPTIC FLOW                                 */
/*                                                                          */
/*                 (Copyright by Joachim Weickert, 8/2014)                  */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 features:
 - explicit scheme
 - Charbonnier diffusivity
*/


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

void alloc_cubix

     (float ****cubix,  /* cubix */
      long  n1,         /* size in direction 1 */
      long  n2,         /* size in direction 2 */
      long  n3)         /* size in direction 3 */

     /* allocates memory for cubix of size n1 * n2 * n3 */


{
long i, j;

*cubix = (float ***) malloc (n1 * sizeof(float **));
if (*cubix == NULL) 
   {
   printf("alloc_cubix: not enough memory available\n");
   exit(1);
   }
for (i=0; i<n1; i++) 
    {
    (*cubix)[i] = (float **) malloc (n2 * sizeof(float *));
    if ((*cubix)[i] == NULL) 
       {
       printf("alloc_cubix: not enough memory available\n");
       exit(1);
       }
    for (j=0; j<n2; j++) 
        {
        (*cubix)[i][j] = (float *) malloc (n3 * sizeof(float));
        if ((*cubix)[i][j] == NULL) 
           {
           printf("alloc_cubix: not enough memory available\n");
           exit(1);
           }
        }
    }
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

/*----------------------------------------------------------------------------*/

void disalloc_cubix

     (float ***cubix,   /* cubix */
      long  n1,         /* size in direction 1 */
      long  n2,         /* size in direction 2 */
      long  n3)         /* size in direction 3 */

     /* disallocates memory for cubix of size n1 * n2 * n3 */

{
long i, j;

for (i=0; i<n1; i++) 
 for (j=0; j<n2; j++) 
     free(cubix[i][j]);

for (i=0; i<n1; i++)
    free(cubix[i]);

free(cubix);

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

/*----------------------------------------------------------------------------*/

float maxf(float a,float b)
{
    return a>b ? a : b;
}

/*---------------------------------------------------------------------------*/

float minf(float a,float b)
{
    return a<b ? a : b;
}

/*---------------------------------------------------------------------------*/

int maxi(int a,int b)
{
  return (a>b) ? a : b;
}

/*---------------------------------------------------------------------------*/

int mini(int a,int b)
{
  return (a<b) ? a : b;
}

/*---------------------------------------------------------------------------*/

int clamped (int min, int val, int max)
{
    assert(min<=max);
    return mini(maxi(min,val),max);
}

/*----------------------------------------------------------------------------*/

int byte_range(int a)
/* restricts number to unsigned char range */
{
    return clamped(0,a,255);
}

/*----------------------------------------------------------------------------*/

void vector_to_RGB

     (float x,   /* x-component */
      float y,   /* y-component */
      int   *R,  /* red component */
      int   *G,  /* green component */
      int   *B)  /* blue component */

/* 
  Computes the color representation of a vector. 
*/

{               
float Pi;          /* pi */
float amp;         /* amplitude (magnitude) */
float phi;         /* phase (angle) */
float alpha, beta; /* weights for linear interpolation */

/* set pi */
Pi = 2.0 * acos(0.0);

/* determine amplitude and phase (cut amp at 1) */
amp = sqrt (x * x + y * y);
if (amp > 1) amp = 1;
if (x == 0.0)
  if (y >= 0.0) phi = 0.5 * Pi;
  else phi = 1.5 * Pi;
else if (x > 0.0)
  if (y >= 0.0) phi = atan (y/x);
  else phi = 2.0 * Pi + atan (y/x);
else phi = Pi + atan (y/x);

phi = phi / 2.0;

// interpolation between red (0) and blue (0.25 * Pi)
if ((phi >= 0.0) && (phi < 0.125 * Pi)) 
   {
   beta  = phi / (0.125 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   *G = (int)floor(amp * (alpha *   0.0 + beta *   0.0));
   *B = (int)floor(amp * (alpha *   0.0 + beta * 255.0));
   }
if ((phi >= 0.125 * Pi) && (phi < 0.25 * Pi)) 
   {
   beta  = (phi-0.125 * Pi) / (0.125 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha * 255.0 + beta *  64.0));
   *G = (int)floor(amp * (alpha *   0.0 + beta *  64.0));
   *B = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   }

// interpolation between blue (0.25 * Pi) and green (0.5 * Pi)
if ((phi >= 0.25 * Pi) && (phi < 0.375 * Pi)) 
   {
   beta  = (phi - 0.25 * Pi) / (0.125 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha *  64.0 + beta *   0.0));
   *G = (int)floor(amp * (alpha *  64.0 + beta * 255.0));
   *B = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   }
if ((phi >= 0.375 * Pi) && (phi < 0.5 * Pi)) 
   {
   beta  = (phi - 0.375 * Pi) / (0.125 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha *   0.0 + beta *   0.0));
   *G = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   *B = (int)floor(amp * (alpha * 255.0 + beta *   0.0));
   }
   
// interpolation between green (0.5 * Pi) and yellow (0.75 * Pi)
if ((phi >= 0.5 * Pi) && (phi < 0.75 * Pi)) 
   {
   beta  = (phi - 0.5 * Pi) / (0.25 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha * 0.0   + beta * 255.0));
   *G = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   *B = (int)floor(amp * (alpha * 0.0   + beta * 0.0));
   }

// interpolation between yellow (0.75 * Pi) and red (Pi)
if ((phi >= 0.75 * Pi) && (phi <= Pi)) 
   {
   beta  = (phi - 0.75 * Pi) / (0.25 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   *G = (int)floor(amp * (alpha * 255.0 + beta *   0.0));
   *B = (int)floor(amp * (alpha * 0.0   + beta *   0.0));
   }

/* check RGB range */
*R = byte_range(*R);
*G = byte_range(*G);
*B = byte_range(*B);

return;

} /* vector_to_RGB */

/*----------------------------------------------------------------------------*/

void flow_to_color

     (float **u,            /* flow field, first channel (input) */
      float **v,            /* flow field, second channel (input) */
      float ***color_img,   /* color representation (output) */
      float max_disp,       /* maximal disparity (set to -1 if not used) */
      int nx,               /* size in x direction */
      int ny)               /* size in y direction */

/* 
  Computes a color representation of a flow field. 
*/

{
int i,j;
int R,G,B;
float maximum_length = 0;

for (i=1; i<=nx; i++) 
 for (j=1; j<=ny; j++) 
     maximum_length = maxf(maximum_length, u[i][j] * u[i][j]
                                         + v[i][j] * v[i][j]);

if(max_disp==-1.0f)
   maximum_length = 1.0f / sqrt(maximum_length);
else 
   maximum_length = 1.0f / max_disp;

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++) 
     {
     if(u[i][j]!=100.0f && v[i][j]!=100.0f) 
        vector_to_RGB(maximum_length*u[i][j],
                      maximum_length*v[i][j],&R,&G,&B);
     else
        R=G=B=128.0f;

     color_img[i][j][0]=(float)R;
     color_img[i][j][1]=(float)G;
     color_img[i][j][2]=(float)B;
     }

return;

} /* flow_to_color */

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

/*----------------------------------------------------------------------------*/

void write_ppm

     (float ***u,          /* colour image, unchanged */
      int   nx,            /* size in x direction */
      int   ny,            /* size in y direction */
      char  *file_name,    /* name of ppm file */
      char  *comments)     /* comment string (set 0 for no comments) */

     /* writes an image into a pgm P5 (greyscale) or ppm P6 (colour) file */

{
FILE           *outimage;  /* output file */
int            i, j, m;    /* loop variables */
float          aux;        /* auxiliary variable */
unsigned char  byte;       /* for data conversion */

/* open file */
outimage = fopen (file_name, "wb");
if (NULL == outimage) 
   {
   printf("Could not open file '%s' for writing, aborting\n", file_name);
   exit(1);
   }

fprintf (outimage, "P6\n");                /* colour format */

if (comments != 0)
   fprintf (outimage, comments);           /* comments */

fprintf (outimage, "%d %d\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");               /* maximal value */

/* write image data */
for (j = 1; j <= ny; j++) 
 for (i = 1; i <= nx; i++) 
  for (m = 0; m < 3; m++) 
      {
      aux = u[i][j][m] + 0.499999;    /* for correct rounding */
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
  
} /* write_ppm */

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

void flow 

     (float    ht,          /* time step size, 0 < ht <= 0.25 */
      long     nx,          /* image dimension in x direction */ 
      long     ny,          /* image dimension in y direction */ 
      float    hx,          /* pixel size in x direction */
      float    hy,          /* pixel size in y direction */
      float    **fx,        /* x derivative of image */
      float    **fy,        /* y derivative of image */
      float    **ft,        /* theta derivative of image */
      float    alpha,       /* smoothness weight */
      float    lambda,      /* contrast parameter */
      float    **u,         /* x component of optic flow */
      float    **v)         /* v component of optic flow */

/* 
 Optic flow iteration by regarding it as a coupled system of
 two nonlinear diffusion-reaction equations.
 Isotropic nonlinear diffusion with explicit discretization.
*/

{
long    i, j;             /* loop variables */
float   rx, ry;           /* mesh ratios, time savers */
float   **u1, **v1;       /* intermediate results */
float   **dc;             /* diffusion coefficient */
float   **rhs;            /* right hand side of linear system */
float   **diag1;          /* summand of diagonal part of system matrix */
float   ux, uy, vx, vy;   /* derivatives */
float   help, rxx, ryy;   /* time savers */
float   two_hx, two_hy;   /* time savers */
float   grad_sqr;         /* |grad(v)|^2 */
float   lambda_sqr;       /* lambda^2 */
float   n, s, w, e, c;    /* stencil coefficients: north, south, west, 
                             east, centre */


/* ---- allocate storage ---- */

alloc_matrix (&u1, nx+2, ny+2);
alloc_matrix (&v1, nx+2, ny+2);
alloc_matrix (&dc, nx+2, ny+2);
alloc_matrix (&rhs, nx+2, ny+2);
alloc_matrix (&diag1, nx+2, ny+2);


/* ---- copy u, v into u1, v1 ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     u1[i][j] = u[i][j];
     v1[i][j] = v[i][j];
     }


/* ---- calculate diffusivity ---- */

two_hx = 2.0 * hx;
two_hy = 2.0 * hy;
lambda_sqr = lambda * lambda;
dummies (u1, nx, ny);
dummies (v1, nx, ny);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* calculate grad(f) */
     ux = (u1[i+1][j] - u1[i-1][j]) / two_hx;
     uy = (u1[i][j+1] - u1[i][j-1]) / two_hy;
     vx = (v1[i+1][j] - v1[i-1][j]) / two_hx;
     vy = (v1[i][j+1] - v1[i][j-1]) / two_hy;
     grad_sqr = ux * ux + uy * uy + vx * vx + vy * vy;
     dc[i][j] = 1.0 / sqrt (1.0 + grad_sqr / lambda_sqr);
     }

/* create reflecting dummy boundaries */
dummies (dc, nx, ny);


/* ---- calculate diffusion of the first equation ---- */

help = ht / alpha;
rxx  = ht / (2.0 * hx * hx);
ryy  = ht / (2.0 * hy * hy);

/* create reflecting dummy boundaries */
dummies (u1, nx, ny);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* SUPPLEMENT CODE */

     /* u[i][j] = ... */
     }


/* ---- calculate diffusion of the second equation ---- */

/* create reflecting dummy boundaries */
dummies (v1, nx, ny);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* SUPPLEMENT CODE */

     /* v[i][j] = ... */
     }


/* ---- disallocate storage ---- */

disalloc_matrix (u1, nx+2, ny+2);
disalloc_matrix (v1, nx+2, ny+2);
disalloc_matrix (dc, nx+2, ny+2);
disalloc_matrix (rhs, nx+2, ny+2);
disalloc_matrix (diag1, nx+2, ny+2);

return;

} /* flow */

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
char   in1[80], in2[80];     /* for reading data */
char   out[80];              /* for reading data */
float  **f1, **f2;           /* images */
float  **fx, **fy, **ft;     /* image derivatives */
float  **u, **v;             /* optic flow components */
float  **w;                  /* optic flow magnitude */
float  ***color_img;         /* color representation of optic flow field */
long   dtype;                /* type of diffusivity */
long   i, j, k;              /* loop variables */
long   kmax;                 /* max. no. of iterations */
long   nx, ny;               /* image size in x, y direction */
float  hx, hy;               /* pixel sizes */
float  ht;                   /* time step size */ 
float  alpha;                /* smoothness weight */
float  lambda;               /* contrast parameter */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("OPTIC FLOW WITH CHARBONNIER SMOOTHNESS TERM, EXPLICIT SCHEME\n\n");
printf ("************************************************************\n\n");
printf ("    Copyright 2014 by Joachim Weickert                      \n");
printf ("    Dept. of Mathematics and Computer Science               \n");
printf ("    Saarland University, Saarbruecken, Germany              \n\n");
printf ("    All rights reserved. Unauthorized usage,                \n");
printf ("    copying, hiring, and selling prohibited.                \n\n");
printf ("    Send bug reports to                                     \n");
printf ("    weickert@mia.uni-saarland.de                            \n\n");
printf ("************************************************************\n\n");


/* ---- read input image f1 (pgm format P5) ---- */

printf ("input image 1 (pgm):                     ");
read_string (in1);
read_pgm_and_allocate_memory (in1, &nx, &ny, &f1);


/* ---- read input image f2 (pgm format P5) ---- */

printf ("input image 2 (pgm):                     ");
read_string (in2);
read_pgm_and_allocate_memory (in2, &nx, &ny, &f2);


/* ---- read parameters ---- */

printf ("smoothnes weight alpha (>0) (float):     ");
read_float (&alpha);

printf ("contrast parameter lambda (>0) (float):  ");
read_float (&lambda);

printf ("time step size (<=0.25) (float):         ");
read_float (&ht);

printf ("number of iterations (>0) (integer):     ");
read_long (&kmax);

printf ("output image (colour coded flow) (ppm):     ");
read_string (out);
printf ("\n");


/* ---- initializations ---- */

/* allocate storage for image derivatives fx, fy, ft */
alloc_matrix (&fx, nx+2, ny+2);
alloc_matrix (&fy, nx+2, ny+2);
alloc_matrix (&ft, nx+2, ny+2);

/* calculate image derivatives fx, fy and ft */
dummies (f1, nx, ny);
dummies (f2, nx, ny);
hx = 1.0;
hy = 1.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     fx[i][j] = (f1[i+1][j] - f1[i-1][j] + f2[i+1][j] - f2[i-1][j]) 
                / (4.0 * hx);
     fy[i][j] = (f1[i][j+1] - f1[i][j-1] + f2[i][j+1] - f2[i][j-1]) 
                / (4.0 * hy);
     ft[i][j] = f2[i][j] - f1[i][j];   /* frame distance 1 assumed */ 
     }
 
/* allocate storage for optic flow and its magnitude */
alloc_matrix (&u, nx+2, ny+2);
alloc_matrix (&v, nx+2, ny+2);
alloc_matrix (&w, nx+2, ny+2);

/* allocate storage for output */
alloc_cubix (&color_img, nx+2, ny+2, 3);

/* initialize (u,v) with 0 */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     u[i][j] = 0.0;
     v[i][j] = 0.0;
     }


/* ---- process image ---- */

for (k=1; k<=kmax; k++)
    {
    /* perform one iteration */
    printf ("iteration number: %5ld \n", k);
    flow (ht, nx, ny, hx, hy, fx, fy, ft, alpha, lambda, u, v);

    /* calculate flow magnitude */
    for (i=1; i<=nx; i++)
     for (j=1; j<=ny; j++)
         w[i][j] = sqrt (u[i][j] * u[i][j] + v[i][j] * v[i][j]);

    /* check maximum, mean, standard deviation of flow magnitude */
    analyse (w, nx, ny, &min, &max, &mean, &std);
    printf ("maximum:       %8.2f \n", max);
    printf ("mean:          %8.2f \n", mean);
    printf ("standard dev.: %8.2f \n\n", std);
    } 


/* ---- write output image (color representation of flow) ---- */

flow_to_color(u, v, color_img, -1, nx, ny);

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# optic flow (color coded), Charbonnier smoothness\n");
comment_line (comments, "# initial image 1:  %s\n", in1);
comment_line (comments, "# initial image 2:  %s\n", in2);
comment_line (comments, "# alpha:            %8.4f\n", alpha);
comment_line (comments, "# lambda:           %8.4f\n", lambda);
comment_line (comments, "# ht:               %8.4f\n", ht);
comment_line (comments, "# iterations:       %8ld\n",  kmax);
comment_line (comments, "# min:              %8.4f\n", min);
comment_line (comments, "# max:              %8.4f\n", max);
comment_line (comments, "# mean:             %8.4f\n", mean);
comment_line (comments, "# standard dev.:    %8.4f\n", std);

/* write image */
write_ppm (color_img, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (f1, nx+2, ny+2);
disalloc_matrix (f2, nx+2, ny+2);
disalloc_matrix (fx, nx+2, ny+2);
disalloc_matrix (fy, nx+2, ny+2);
disalloc_matrix (ft, nx+2, ny+2);
disalloc_matrix (u,  nx+2, ny+2);
disalloc_matrix (v,  nx+2, ny+2);
disalloc_matrix (w,  nx+2, ny+2);
disalloc_cubix (color_img, nx+2, ny+2, 3);

return(0);
}

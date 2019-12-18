#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define float double

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                            OSMOSIS FILTERING                             */
/*                                                                          */
/*                  (Copyright Joachim Weickert, 6/2013)                    */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 features:
 - for greyscale images
 - explicit scheme
 - test, if an image can be reconstructed by its canonical drift vectors 
*/


/*--------------------------------------------------------------------------*/

void alloc_vector

     (float **vector,   /* vector */
      long  n)          /* size */

     /* allocates storage for a vector of size n */


{
*vector = (float *) malloc (n * sizeof(float));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough storage available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (float ***matrix,  /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* allocates storage for matrix of size nx * ny */


{
long i;

*matrix = (float **) malloc (nx * sizeof(float *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough storage available\n");
   exit(1);
   }
for (i=0; i<nx; i++)
    {
    (*matrix)[i] = (float *) malloc (ny * sizeof(float));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_matrix: not enough storage available\n");
       exit(1);
       }
    }
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

     (float *vector,    /* vector */
      long  n)          /* size */

     /* disallocates storage for a vector of size n */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (float **matrix,   /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* disallocates storage for matrix of size nx * ny */

{
long i;
for (i=0; i<nx; i++)
    free(matrix[i]);
free(matrix);
return;
}

/*--------------------------------------------------------------------------*/

void dummies

     (float **v,        /* image matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

/* creates dummy boundaries by mirroring */

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    v[i][0]    = v[i][1];
    v[i][ny+1] = v[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    v[0][j]    = v[1][j];
    v[nx+1][j] = v[nx][j];
    }
return;
}

/* ---------------------------------------------------------------------- */

void canonical_drift_vectors 

     (float    **v,     /* guidance image, unchanged */
      long     nx,      /* image dimension in x direction */ 
      long     ny,      /* image dimension in y direction */ 
      float    hx,      /* pixel size in x direction */
      float    hy,      /* pixel size in y direction */
      float    **d1,    /* drift vector, x component in [i+1/2,j], output */
      float    **d2)    /* drift vector, y component in [i,j+1/2], output */

/*
 computes the canonical drift vector field that allows to reconstruct the 
 guidance image up to a multiplicative constant
*/

{
long    i, j;             /* loop variables */


/* ---- dummy boundaries for v ---- */

dummies (v, nx, ny);


/* ---- initialise drift vector field with 0 ---- */

for (i=0; i<=nx+1; i++)
 for (j=0; j<=ny+1; j++)
     d1[i][j] = d2[i][j] = 0.0;


/* ---- compute x component of canonical drift vector field ---- */

/* index [i,j] refers to intergrid location [i+1/2,j] */
    
/* SUPPLEMENT CODE */


/* ---- compute y component of canonical drift vector field ---- */

/* index [i,j] refers to intergrid location [i,j+1/2] */
    
/* SUPPLEMENT CODE */


/* ---- modification at the shadow boundaries between i=128 and i=129 ---- */

/* SUPPLEMENT CODE */

return;

} /* canonical_drift_vectors */

/* ------------------------------------------------------------------------ */

void osmosis_weights 

     (float    ht,      /* time step size, 0 < ht < 0.125 */
      long     nx,      /* image dimension in x direction */
      long     ny,      /* image dimension in y direction */
      float    hx,      /* pixel size in x direction */
      float    hy,      /* pixel size in y direction */
      float    **d1,    /* drift vector, x component in [i+1/2,j], unchanged */
      float    **d2,    /* drift vector, y component in [i,j+1/2], unchanged */
      float    **woo,   /* osmosis weight for pixel [i][j],   output */
      float    **wpo,   /* osmosis weight for pixel [i+1][j], output */
      float    **wmo,   /* osmosis weight for pixel [i-1][j], output */
      float    **wop,   /* osmosis weight for pixel [i][j+1], output */
      float    **wom)   /* osmosis weight for pixel [i][j-1], output */

/*
 computes the weights for osmosis filtering
*/

{
long    i, j;             /* loop variables */
float   rx, rxx;          /* time savers */
float   ry, ryy;          /* time savers */


/* ---- initialise all osmosis weights ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     woo[i][j] = 1.0;
     wpo[i][j] = wmo[i][j] = wop[i][j] = wom[i][j] = 0.0; 
     }


/* ---- specify them from the drift vector field ---- */

/* compute time savers */
rx  = ht / (2.0 * hx);
ry  = ht / (2.0 * hy);
rxx = ht / (hx * hx);
ryy = ht / (hy * hy);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* osmosis weight for pixel [i][j] */
     woo[i][j] = 1.0 - 2.0 * (rxx + ryy)
                 + rx * (d1[i-1][j] - d1[i][j])
                 + ry * (d2[i][j-1] - d2[i][j]);

     /* osmosis weight for pixel [i+1][j] */
     wpo[i][j] = rxx - rx * d1[i][j];

     /* osmosis weight for pixel [i-1][j] */
     wmo[i][j] = rxx + rx * d1[i-1][j];

     /* osmosis weight for pixel [i][j+1] */
     wop[i][j] = ryy - ry * d2[i][j];

     /* osmosis weight for pixel [i][j-1] */
     wom[i][j] = ryy + ry * d2[i][j-1];
     }

return;  /* weights */

}

/* ------------------------------------------------------------------------ */

void osmosis 

     (long     nx,      /* image dimension in x direction */ 
      long     ny,      /* image dimension in y direction */ 
      float    **woo,   /* osmosis weight for pixel [i][j],   output */
      float    **wpo,   /* osmosis weight for pixel [i+1][j], output */
      float    **wmo,   /* osmosis weight for pixel [i-1][j], output */
      float    **wop,   /* osmosis weight for pixel [i][j+1], output */
      float    **wom,   /* osmosis weight for pixel [i][j-1], output */
      float    **u)     /* input: original image;  output: filtered */

/* 
 Osmosis scheme. Explicit discretisation.
*/

{
long    i, j;             /* loop variables */
float   **f;              /* work copy of u */
      

/* ---- allocate storage ---- */

alloc_matrix (&f, nx+2, ny+2);


/* ---- copy u into f ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];


/* ---- dummy boundaries for f ---- */

dummies (f, nx, ny);


/* ---- calculate explicit osmosis of u ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = woo[i][j] * f[i][j] + 
               wpo[i][j] * f[i+1][j] +
               wmo[i][j] * f[i-1][j] +
               wop[i][j] * f[i][j+1] +
               wom[i][j] * f[i][j-1];


/* ---- disallocate storage ---- */

disalloc_matrix (f, nx+2, ny+2);

return;

} /* osmosis */

/* ---------------------------------------------------------------------- */

void analyse

     (float   **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in x direction */
      float   *min,        /* minimum, output */
      float   *max,        /* maximum, output */
      float   *mean,       /* mean, output */
      float   *std)        /* standard deviation, output */

/*
 calculates minimum, maximum, mean and standard deviation of an image u
*/

{
long    i, j;       /* loop variables */
float   help;       /* auxiliary variable */
double  help2;      /* auxiliary variable */

*min  = u[1][1];
*max  = u[1][1];
help2 = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     if (u[i][j] < *min) *min = u[i][j];
     if (u[i][j] > *max) *max = u[i][j];
     help2 = help2 + (double)u[i][j];
     }
*mean = (float)help2 / (nx * ny);

*std = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     help  = u[i][j] - *mean;
     *std = *std + help * help;
     }
*std = sqrt (*std / (nx * ny));

return;

} /* analyse */

/*--------------------------------------------------------------------------*/

int main ()

{
char   row[80];              /* for reading data */
char   in1[80], in2[80];     /* for reading data */
char   out[80];              /* for reading data */
float  **u;                  /* evolving image */
float  **v;                  /* guidance image */
float  **d1;                 /* drift vector field, x component */
float  **d2;                 /* drift vector field, y component */
float  **woo;                /* osmosis weight for pixel [i][j] */
float  **wpo;                /* osmosis weight for pixel [i+1][j] */
float  **wmo;                /* osmosis weight for pixel [i-1][j] */
float  **wop;                /* osmosis weight for pixel [i][j+1] */
float  **wom;                /* osmosis weight for pixel [i][j-1] */
long   i, j, k;              /* loop variables */
long   nx, ny;               /* image size in x, y direction */
FILE   *inimage, *outimage;  /* input file, output file */
float  ht;                   /* time step size */
float  offset;               /* greyscale offset */
long   kmax;                 /* largest iteration number */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
float  aux;                  /* auxiliary variable */
unsigned char byte;          /* for data conversion */


printf("\n");
printf("OSMOSIS FILTERING, EXPLICIT SCHEME\n");
printf("RECONSTRUCTION FROM THE CANONICAL DRIFT VECTOR FIELD\n\n");
printf("*************************************************\n\n");
printf("    Copyright 2013 by Joachim Weickert           \n");
printf("    Dept. of Mathematics and Computer Science    \n");
printf("    Saarland University, Germany                 \n\n");
printf("    All rights reserved. Unauthorized usage,     \n");
printf("    copying, hiring, and selling prohibited.     \n\n");
printf("    Send bug reports to                          \n");
printf("    weickert@mia.uni-saarland.de                 \n\n");
printf("*************************************************\n\n");


/* ---- read initial image (pgm format P5) ---- */

/* read image name */
printf("initial image f:                  ");
gets (in1);

/* open pgm file and read header */
inimage = fopen(in1,"r");
fgets (row, 300, inimage);
fgets (row, 300, inimage);
while (row[0]=='#') fgets(row, 300, inimage);
sscanf (row, "%ld %ld", &nx, &ny);
fgets (row, 300, inimage);

/* allocate storage */
alloc_matrix (&u, nx+2, ny+2);

/* read image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     u[i][j] = (float) getc (inimage);
fclose(inimage);


/* ---- read guidance image (pgm format P5) ---- */

/* read image name */
printf("guidance image v:                 ");
gets (in2);

/* open pgm file and read header */
inimage = fopen(in2,"r");
fgets (row, 300, inimage);
fgets (row, 300, inimage);
while (row[0]=='#') fgets(row, 300, inimage);
sscanf (row, "%ld %ld", &nx, &ny);
fgets (row, 300, inimage);

/* allocate storage */
alloc_matrix (&v, nx+2, ny+2);

/* read image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     v[i][j] = (float) getc (inimage);
fclose(inimage);


/* ---- read other parameters ---- */

printf("time step size (<0.125):          ");
gets(row);  sscanf(row, "%lf", &ht);
printf("number of iterations:             ");
gets(row);  sscanf(row, "%ld", &kmax);
printf("greyscale offset (>0.0):          ");
gets(row);  sscanf(row, "%lf", &offset);
printf("output image:                     ");
gets(out);
printf("\n");


/* ---- allocate storage ---- */

alloc_matrix (&d1,  nx+2, ny+2);
alloc_matrix (&d2,  nx+2, ny+2);
alloc_matrix (&woo, nx+2, ny+2);
alloc_matrix (&wpo, nx+2, ny+2);
alloc_matrix (&wmo, nx+2, ny+2);
alloc_matrix (&wop, nx+2, ny+2);
alloc_matrix (&wom, nx+2, ny+2);


/* ---- process image ---- */

/* add offset in order to make data positive */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     u[i][j] = u[i][j] + offset;
     v[i][j] = v[i][j] + offset;
     }

/* analyse initial image */
printf("initial image\n");
analyse (u, nx, ny, &min, &max, &mean, &std);
printf("minimum:       %10.2lf \n", min);
printf("maximum:       %10.2lf \n", max);
printf("mean:          %10.2lf \n", mean);
printf("std. dev.:     %10.2lf \n\n", std);

/* compute canonical drift vectors of the guidance image */
canonical_drift_vectors (v, nx, ny, 1.0, 1.0, d1, d2);

/* compute resulting osmosis weights */
osmosis_weights (ht, nx, ny, 1.0, 1.0, d1, d2, woo, wpo, wmo, wop, wom);

/* perform kmax osmosis iterations */
for (k=1; k<=kmax; k++)
    {
    /* perform one iteration */
    osmosis (nx, ny, woo, wpo, wmo, wop, wom, u);

    /* check minimum, maximum, mean, standard deviation */
    analyse (u, nx, ny, &min, &max, &mean, &std);
    printf("iteration number: %7ld \n", k);
    printf("minimum:       %10.2lf \n", min);
    printf("maximum:       %10.2lf \n", max);
    printf("mean:          %10.2lf \n", mean);
    printf("std. dev.:     %10.2lf \n\n", std);
    } /* for */

/* subtract offset */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = u[i][j] - offset;


/* ---- write output image (pgm format P5) ---- */

/* open file and write header (incl. filter parameters) */
outimage = fopen (out, "w");
fprintf (outimage, "P5 \n");
fprintf (outimage, "# linear osmosis, explicit scheme\n"); 
fprintf (outimage, "# initial image:  %s\n", in1);
fprintf (outimage, "# guidance image: %s\n", in2);
fprintf (outimage, "# ht:             %8.4lf\n", ht);
fprintf (outimage, "# iterations:     %8ld\n",   kmax);
fprintf (outimage, "# offset:         %8.2lf\n", offset);
fprintf (outimage, "# minimum:        %8.2lf\n", min);
fprintf (outimage, "# maximum:        %8.2lf\n", max);
fprintf (outimage, "# mean:           %8.2lf\n", mean);
fprintf (outimage, "# std. dev.:      %8.2lf\n", std);
fprintf (outimage, "%ld %ld \n255\n", nx, ny);

/* write image data and close file */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     aux = u[i][j] + 0.49999;
     if (aux < 0.0)
        byte = (unsigned char)(0.0);
     else if (aux > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(aux);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }
fclose(outimage);
printf("output image %s successfully written\n\n", out);


/* ---- disallocate storage ---- */

disalloc_matrix (u,   nx+2, ny+2);
disalloc_matrix (v,   nx+2, ny+2);
disalloc_matrix (d1,  nx+2, ny+2);
disalloc_matrix (d2,  nx+2, ny+2);
disalloc_matrix (woo, nx+2, ny+2);
disalloc_matrix (wpo, nx+2, ny+2);
disalloc_matrix (wmo, nx+2, ny+2);
disalloc_matrix (wop, nx+2, ny+2);
disalloc_matrix (wom, nx+2, ny+2);

return(0);
}

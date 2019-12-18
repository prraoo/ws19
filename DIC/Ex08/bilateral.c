#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                       ITERATED BILATERAL FILTERING                       */
/*                                                                          */
/*   (Copyright by Markus Mainberger, 5/2009 and Joachim Weickert, 8/2014)  */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 iterated bilateral filtering
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

void bilateral 

     (long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    theta_t,   /* standard deviation - tonal weight */
      float    theta_s,   /* standard deviation - spatial weight */
      float    **u)       /* input: original image; output: processed */

/*
 bilateral filtering
*/

{
float h_tonal = -1.0f / (theta_t * theta_t);    /* time saver */
float h_spatial = -1.0f / (theta_s * theta_s);  /* time saver */
float spatial_precision = 3.0f * theta_s;       /* time saver */

float  **u_old;               /* evolving image, prev. iteration */
long   i1, i2, j1, j2, k;     /* loop variables */
long   l1, l2, r1, r2;        /* loop variables for domain restriction */
float  g, w;                  /* auxiliary variables for computing
                                tonal/ spatial weights */
float numerator, denominator; /* auxiliary variables for determining the 
                                 numerator/ denominator of the solution
                                 in pixel (i1,i2) */


/* ---- allocate storage for u_old ---- */

alloc_matrix (&u_old,  nx+2, ny+2);

/* copy image from previous iteration */
for (i2=1; i2<=ny; i2++)
 for (i1=1; i1<=nx; i1++)
     u_old[i1][i2] = u[i1][i2];


/* ---- process image ---- */

/* compute the solution for all image pixels (i1,i2) */
for (i1=1; i1<=nx; i1++)
 for (i2=1; i2<=ny; i2++)
     {
     numerator = 0.0f;
     denominator = 0.0f;

     /* restrict domain to save compuation time */
     l1 = (i1-spatial_precision);
     l1 = (l1 < 1) ? 1 : l1;

     r1 = (i1+spatial_precision);
     r1 = (r1 > nx) ? nx : r1;

     l2 = (i2-spatial_precision);
     l2 = (l2 < 1) ? 1 : l2;

     r2 = (i2+spatial_precision);
     r2 = (r2 > ny) ? ny : r2;

     /* sum over all j in J */
     for (j1 = l1; j1 <= r1; j1++)
      for (j2 = l2; j2 <= r2; j2++)
          {
          /* compute tonal weight */
          g = expf((u_old[i1][i2]-u_old[j1][j2]) * (u_old[i1][i2]-u_old[j1][j2]) * h_tonal);


          /* compute spatial weight */
          w = expf((((i1-j1)*(i1-j1)) + ((i2-j2)*(i2-j2))) * h_spatial);

          /* add weight of corresponding pixel (j1,j2) to numerator
             and denominator */
          numerator +=  g*w*u_old[j1][j2];
          denominator +=  g*w;
          }

          /* compute solution in pixel (i1,i2) */
          u[i1][i2] = numerator / denominator;
     }


/* ---- disallocate storage for u_old ---- */

disalloc_matrix (u_old,  nx+2, ny+2);

return;

} /* bilateral */

/*--------------------------------------------------------------------------*/

void analyse

     (float   **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in x direction */
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
long   k;                    /* loop variable */
long   nx, ny;               /* image size in x, y direction */ 
long   kmax;                 /* number of iterations */
float  theta_t;              /* standard deviation - tonal weight */
float  theta_s;              /* standard deviation - spatial weight */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("ITERATED BILATERAL FILTERING\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2009 by Markus Mainberger           \n");
printf ("    and 2014 by Joachim Weickert                  \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorized usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                        ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &u);


/* ---- read parameters ---- */

printf ("tonal parameter theta_t (>0) (float):     ");
read_float (&theta_t);

printf ("spatial parameter theta_s (>0) (float):   ");
read_float (&theta_s);

printf ("number of iterations (>0) (integer):      ");
read_long (&kmax);

printf ("output image (pgm):                       ");
read_string (out);
printf ("\n");


/* ---- process image ---- */

for (k=1; k<=kmax; k++)
    {
    /* perform one iteration */
    printf ("iteration number: %5ld \n\n", k);
    bilateral (nx, ny, theta_t, theta_s, u);
    
    
    /* ---- analyse filtered image ---- */

    /* check minimum, maximum, mean, standard deviation */
    analyse (u, nx, ny, &min, &max, &mean, &std);
    printf ("minimum:       %8.2f \n", min);
    printf ("maximum:       %8.2f \n", max);
    printf ("mean:          %8.2f \n", mean);
    printf ("standard dev.: %8.2f \n\n", std);
    }


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# iterated bilateral filtering\n");
comment_line (comments, "# theta_s:       %8.4f\n", theta_s);
comment_line (comments, "# theta_t:       %8.4f\n", theta_t);
comment_line (comments, "# iterations:    %8ld\n", kmax);
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

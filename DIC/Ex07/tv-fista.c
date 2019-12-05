#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*             					FISTA FOR TV REGULARISATION             						*/
/*                                                                          */
/*          (Copyright by Sven Grewenig, Simon Setzer, 5/2013 and           */
/*                         Joachim Weickert, 8/2014)                        */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 TV regularisation:
 - forward-backward splitting
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

void dummies_Neumann
 
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

void proj

     (float   alpha,       /* regularisation weight */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in y direction */
      float   **dx,        /* image derivatives in x dir. */
      float   **dy)        /* image derivatives in y dir. */

/* computes the projection into the set C */

{
long    i, j;       /* loop variables */ 
float   magn;       /* gradient magnitude */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     magn = sqrtf (dx[i][j] * dx[i][j] + dy[i][j] * dy[i][j]);

     /* compute the projection if the gradient magnitude is
        larger than alpha */
     if (magn > alpha)
        {
        dx[i][j] = dx[i][j] * alpha / magn;
        dy[i][j] = dy[i][j] * alpha / magn;
        }
     }

return;

}  /* proj */

/*--------------------------------------------------------------------------*/

void mult_D

     (long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in y direction */
      float   **u,         /* input vector */
      float   **dx,        /* matrix-vector product, x direction */
      float   **dy)        /* matrix-vector product, y direction */

/* 
  computes the matrix-vector product D*u, i.e. computes the forward 
  differences w.r.t. each direction 
*/

{
long    i, j;       /* loop variables */

/* create Neumann dummy boundaries */
dummies_Neumann (u, nx, ny);

/* forward differences */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     dx[i][j] = u[i+1][j  ] - u[i][j];
     dy[i][j] = u[i  ][j+1] - u[i][j];
     }

return;

}  /* mult_D */

/*--------------------------------------------------------------------------*/

void mult_DT

     (long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in y direction */
      float   **u,         /* resulting vector */
      float   **dx,        /* derivative vector, x direction */
      float   **dy)        /* derivative vector, y direction */

/* 
  computes the matrix-vector product u = D^T*d, 
  i.e. computes the negative backward differences 
*/

{
long    i, j;               /* loop variables */
float   dx_aux, dy_aux;     /* auxiliary variables */

for (i=1; i<=nx; i++)
	for (j=1; j<=ny; j++){
		if(i ==1){
			dx_aux = - dx[i][j];
		}else if (i == nx) {
			dx_aux = + dx[i-1][j];
		}else{
			dx_aux = dx[i-1][j] - dx[i][j];
		}
		if(j == 1){
			dy_aux = -dy[i][j];
		}else if (j == ny){
			dy_aux = dy[i][j-1];
		}else{
			dy_aux = dy[i][j-1] - dy[i][j];
		}
		u[i][j] = dx_aux + dy_aux;
	}
return;

}  /* mult_DT */

/*--------------------------------------------------------------------------*/

void FISTA

     (long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      long     kmax,      /* number of iterations */
      float    alpha,     /* regularisation weight */
      float    tau,       /* time step size */
      float    **u)       /* input: original image; output: processed */

/*
 fast iterative shrinkage thresholding algorithm (FISTA)
*/

{
long    i, j, k;    /* loop variables */
float   **b_oldx;   /* old iteration b^k (x) */
float   **b_oldy;   /* old iteration b^k (y) */
float   **b_newx;   /* new iteration b^{k+1} (x) */
float   **b_newy;   /* new iteration b^{k+1} (y) */
float   **bp_oldx;  /* old iteration b^k (x), proj. */
float   **bp_oldy;  /* old iteration b^k (y), proj. */
float   **bp_newx;  /* new iteration b^{k+1} (x), proj. */
float   **bp_newy;  /* new iteration b^{k+1} (y), proj. */
float   **v;        /* auxiliary variable */
float   theta;      /* step size */
float   tk, tk1;    /* step sizes */


/* ---- allocate storage ---- */

alloc_matrix (&b_oldx,  nx+2, ny+2);
alloc_matrix (&b_oldy,  nx+2, ny+2);
alloc_matrix (&b_newx,  nx+2, ny+2);
alloc_matrix (&b_newy,  nx+2, ny+2);
alloc_matrix (&bp_oldx,  nx+2, ny+2);
alloc_matrix (&bp_oldy,  nx+2, ny+2);
alloc_matrix (&bp_newx,  nx+2, ny+2);
alloc_matrix (&bp_newy,  nx+2, ny+2);
alloc_matrix (&v,       nx+2, ny+2);


/* ---- initialisations ---- */

/* compute b^0 = D*u */
mult_D (nx, ny, u, b_newx, b_newy);

/* t_{-1} = 0, t_0 = 1 */
tk  = 0.0;
tk1 = 1.0;

/* b^{k+1} = 0 */
for (i=0; i<nx+2; i++)
 for (j=0; j<ny+2; j++)
     bp_newx[i][j] = bp_newy[i][j] = 0.0;


/* ---- process image ---- */

for (k=0; k<kmax; k++)
{
	
	/* update step sizes tau and theta */

	/*
		SUPPLEMENT CODE
	*/
      tk=tk1;
      tk1 = (1 + sqrt(1 + 4*tk*tk))/2;
      theta = (tk-1)/tk1;

	/* create b^{k} and projected b^{k} */
	for (i=1; i<=nx; i++)
		for (j=1; j<=ny; j++)
		{
			/*
			SUPPLEMENT CODE
			*/
            b_oldx[i][j] = b_newx[i][j];
            b_oldy[i][j] = b_newy[i][j];
            
            bp_oldx[i][j] = bp_newx[i][j];
            bp_oldy[i][j] = bp_newy[i][j];
		}

	/* compute v = D^T*b^k */
	
	/*
	 SUPPLEMENT CODE
	 */
	
      mult_DT (nx, ny, v, b_oldx, b_oldy);

	/* v = D^T*b^k - u */
	for (i=1; i<=nx; i++)
		for (j=1; j<=ny; j++)
			v[i][j] = v[i][j] - u[i][j];

	/* compute b= D*v^k */
	
	/*
	 SUPPLEMENT CODE
	 */

      mult_D (nx, ny, v, b_newx, b_newy);
      
	/* compute the vector to be projected */
	for (i=1; i<=nx; i++)
		for (j=1; j<=ny; j++)
		{
			/*
				SUPPLEMENT CODE
			*/
            bp_newx[i][j]  = b_oldx[i][j] - tau * b_newx[i][j];
            bp_newy[i][j]  = b_oldy[i][j] - tau * b_newy[i][j];
		}

	/* get projection b^{k+1} */
	
	/*
	 SUPPLEMENT CODE
	 */
      proj (alpha, nx, ny, bp_newx, bp_newy);

	/* relaxation step */
	for (i=1; i<nx+1; i++)
		for (j=1; j<ny+1; j++)
		{
			/*
			SUPPLEMENT CODE
			*/
            
            b_newx[i][j]  = bp_newx[i][j] + theta*(bp_newx[i][j] - bp_oldx[i][j]);
            b_newy[i][j]  = bp_newy[i][j] + theta*(bp_newy[i][j] - bp_oldy[i][j]);
		}
} /* for k */

	/*compute the solution  u - D^T*b ---- */

	mult_DT (nx, ny, v, b_newx, b_newy);

	for (i=1; i<=nx; i++)
		for (j=1; j<=ny; j++)
			u[i][j] = u[i][j] - v[i][j];


/* ---- disallocate storage ---- */

disalloc_matrix (b_oldx,  nx+2, ny+2);
disalloc_matrix (b_oldy,  nx+2, ny+2);
disalloc_matrix (b_newx,  nx+2, ny+2);
disalloc_matrix (b_newy,  nx+2, ny+2);
disalloc_matrix (bp_oldx,  nx+2, ny+2);
disalloc_matrix (bp_oldy,  nx+2, ny+2);
disalloc_matrix (bp_newx,  nx+2, ny+2);
disalloc_matrix (bp_newy,  nx+2, ny+2);
disalloc_matrix (v,       nx+2, ny+2);

return;

} /* FISTA */

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
long   nx, ny;               /* image size in x, y direction */ 
long   kmax;                 /* number of iterations */
float  alpha;                /* regularisation parameter */
float  tau;                  /* time step size */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("TV REGULARISATION WITH FISTA ALGORITHM            \n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2013 by Sven Grewenig, Simon Setzer \n");
printf ("    and 2014 by Joachim Weickert                  \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorized usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                     ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &u);


/* ---- read parameters ---- */

printf ("regularisation parameter (>0) (float): ");
read_float (&alpha);

printf ("time step size (float):                ");
read_float (&tau);

printf ("number of iterations (integer):        ");
read_long (&kmax);

printf ("output image (pgm):                    ");
read_string (out);
printf ("\n");


/* ---- process image ---- */

FISTA (nx, ny, kmax, alpha, tau, u);


/* ---- analyse filtered image ---- */

analyse (u, nx, ny, &min, &max, &mean, &std);
printf ("filtered image:\n");
printf ("minimum:       %8.2f \n", min);
printf ("maximum:       %8.2f \n", max);
printf ("mean:          %8.2f \n", mean);
printf ("standard dev.: %8.2f \n\n", std);


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# TV regularisation with FISTA\n");
comment_line (comments, "# alpha:         %8.4f\n", alpha);
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                 TV RESTORATION WITH THE KACANOV METHOD                   */
/*                                                                          */
/*                       (Joachim Weickert, 5/2013)                         */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

void alloc_vector

     (double **vector,   /* vector */
      long  n)          /* size */

     /* allocates storage for a vector of size n */


{
*vector = (double *) malloc (n * sizeof(double));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough storage available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (double ***matrix,  /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* allocates storage for matrix of size nx * ny */


{
long i;

*matrix = (double **) malloc (nx * sizeof(double *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough storage available\n");
   exit(1);
   }
for (i=0; i<nx; i++)
    {
    (*matrix)[i] = (double *) malloc (ny * sizeof(double));
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

     (double *vector,    /* vector */
      long  n)          /* size */

     /* disallocates storage for a vector of size n */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (double **matrix,   /* matrix */
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

/* ----------------------------------------------------------------------- */

void read_string

(char *v)         /* string to be read */

/*
 r eads a long value v                                                         *
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
 r eads a long value v                                                         *
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

void read_double

(double *v)         /* value to be read */

/*
 r eads a double value v                                                        *
 */

{
	char   row[80];    /* string for reading data */
	
	fgets (row, 80, stdin);
	if (row[strlen(row)-1] == '\n')
		row[strlen(row)-1] = 0;
	sscanf(row, "%lf", &*v);
	return;
}

/*--------------------------------------------------------------------------*/

void read_pgm_and_allocate_memory

(const char  *file_name,    /* name of pgm file */ 
 long        *nx,           /* image size in x direction, output */
 long        *ny,           /* image size in y direction, output */
 double       ***u)          /* image, output */   

/* 
 r eads a greyscale image that has been encoded in pgm format P5;             *
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
			(*u)[i][j] = (double) getc(inimage);
		
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
 A dd a line to the comment string comment. The string line can contain plain *
 text and format characters that are compatible with sprintf.
 Example call: print_comment_line(comment,"Text %f %d",double_var,int_var);
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

(double  **u,          /* image, unchanged */ 
 long   nx,           /* image size in x direction */
 long   ny,           /* image size in y direction */
 char   *file_name,   /* name of pgm file */
 char   *comments)    /* comment string (set 0 for no comments) */

/* 
 w rites a greyscale image into a pgm P5 file;                                *
 */

{
	FILE           *outimage;  /* output file */
	long           i, j;       /* loop variables */
	double          aux;        /* auxiliary variable */
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

     (double **v,        /* image matrix */
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

/* ----------------------------------------------------------------------- */

void neg_dummies

     (double **v,        /* image matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

/* creates negative dummy boundaries by mirroring */

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    v[i][0]    = - v[i][1];
    v[i][ny+1] = - v[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    v[0][j]    = - v[1][j];
    v[nx+1][j] = - v[nx][j];
    }
return;
}

/* ---------------------------------------------------------------------- */

void diffusivity 

     (long     nx,        /* image dimension in x direction */ 
      long     ny,        /* image dimension in y direction */ 
      double    hx,        /* pixel size in x direction */
      double    hy,        /* pixel size in y direction */
      double    eps,       /* parameter epsilon */
      double    **u,       /* original image, unchanged */
      double    **dc)      /* diffusivity, output */


/* 
 Computes epsilon-regularised TV diffusivities. 
*/

{
long    i, j;                 /* loop variables */
double   gx, gy;               /* derivatives */
double   sqr_grad;             /* |grad(u)|^2 */
double   sqr_eps;              /* epsilon * epsilon */
double   aux1, aux2;           /* time savers */

dummies (u, nx, ny);

aux1 = 1.0 / (2.0 * hx);
aux2 = 1.0 / (2.0 * hy);
sqr_eps = eps * eps;

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     gx = aux1 * (u[i+1][j] - u[i-1][j]);
     gy = aux2 * (u[i][j+1] - u[i][j-1]);
     sqr_grad = gx * gx + gy * gy;
     dc[i][j] = 1.0 / sqrt (sqr_eps + sqr_grad);
     }

return;

} /* diffusivity */

/* ---------------------------------------------------------------------- */

void gauss_seidel

    (double   alpha,  /* regularisation parameter */
     long    nx,     /* image dimension in x direction */ 
     long    ny,     /* image dimension in y direction */ 
     double   hx,     /* pixel size in x direction */
     double   hy,     /* pixel size in y direction */
     double   **g,    /* old image, unchanged on exit */
     double   **dc,   /* diffusivity, with boundary mod's on exit */
     double   **f)    /* old and new solution */


/*
 Creates one Gauss-Seidel iteration.
*/

{
long   i, j;                  /* loop variables */
double  wxp, wyp, wxm, wym;    /* weights */
double  r1, r2;                /* time saver */

dummies (f, nx, ny);
neg_dummies (dc, nx, ny);

r1 = alpha / (2.0 * hx);   
r2 = alpha / (2.0 * hy);   

for (i=1; i<=nx; i++) 
 for (j=1; j<=ny; j++) 
     { 
     wxp = dc[i+1][j] + dc[i][j]; 
     wxm = dc[i-1][j] + dc[i][j]; 
     wyp = dc[i][j+1] + dc[i][j]; 
     wym = dc[i][j-1] + dc[i][j]; 

     f[i][j] = (g[i][j] + r1 * (wxp * f[i+1][j] + wxm * f[i-1][j])
                        + r2 * (wyp * f[i][j+1] + wym * f[i][j-1]))
               / (1.0 + r1 * (wxp + wxm) + r2 * (wyp + wym)); 
     } /* for i,j */ 

return;

} /* gauss_seidel */

/*--------------------------------------------------------------------------*/

void analyse

     (double   **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in x direction */
      double   *min,        /* minimum, output */
      double   *max,        /* maximum, output */
      double   *mean,       /* mean, output */
      double   *vari)       /* variance, output */

/*
 calculates minimum, maximum, mean and variance of an image u
*/

{
long    i, j;       /* loop variables */
double   help;       /* auxiliary variable */
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
*mean = (double)help2 / (nx * ny);

*vari = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     help  = u[i][j] - *mean;
     *vari = *vari + help * help;
     }
*vari = *vari / (nx * ny);

return;

} /* analyse */

/*--------------------------------------------------------------------------*/

int main ()

{
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
double  **u;                  /* evolving image */
double  **dc;                 /* diffusivity */
double  **g;                  /* original image */
long   i, j, k, m;           /* loop variables */
long   nx, ny;               /* image size in x, y direction */
double  hx, hy;               /* pixel size in x, y direction */
double  alpha;                /* regularisation parameter */
double  eps;                  /* parameter epsilon */
long   maxfp;                /* number of fixed point iterations */
long   maxgs;                /* mumber of Gauss-Seidel iterations */
double  max, min;             /* largest, smallest grey value */
double  mean;                 /* average grey value */
double  vari;                 /* variance */
char   comments[1600];			 /* comments */

printf("\n");
printf("\n");
printf("TV IMAGE RESTORATION WITH THE KACANOV METHOD\n\n");
printf("*********************************************************\n\n");
printf("    Copyright 2013 by Joachim Weickert                   \n");
printf("    Faculty of Mathematics and Computer Science          \n");
printf("    Saarland University, Saarbruecken, Germany           \n\n");
printf("    All rights reserved. Unauthorized usage, copying,    \n");
printf("    hiring, and selling prohibited.                      \n\n");
printf("    Send bug reports to                                  \n");
printf("    weickert@mia.uni-saarland.de                         \n\n");
printf("*********************************************************\n\n");
printf("- outer iterations:  fixed point\n");
printf("- inner iterations:  Gauss-Seidel\n\n");


/* ---- read input image (pgm format P5) ---- */

/* read image name */
printf("input image:                         ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &u);



/* ---- read other parameters ---- */

printf("regul. parameter alpha (>=0):        ");
read_double (&alpha);

printf("parameter epsilon (>0):              ");
read_double (&eps);

printf("number of fixed point iterations:    ");
read_long (&maxfp);

printf("number of Gauss-Seidel iterations:   ");
read_long (&maxgs);

printf("output image:                        ");
read_string (out);
printf("\n");

/* allocate storage */
alloc_matrix (&dc, nx+2, ny+2);
alloc_matrix (&g,  nx+2, ny+2);


/* ---- analyse initial image ---- */

hx = 1.0;
hy = 1.0;
analyse (u, nx, ny, &min, &max, &mean, &vari);

printf("initial image \n");
printf("minimum:                 %8.2lf \n", min);
printf("maximum:                 %8.2lf \n", max);
printf("mean:                    %8.2lf \n", mean);
printf("variance:                %8.2lf \n\n", vari);


/* ---- process image ---- */

for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     g[i][j] = u[i][j];

for (m=1; m<=maxfp; m++)
  
    /* ---- outer iterations: fixed point ---- */ 

    {     
    printf("fixed point iteration:   %8ld\n", m);

    /* update diffusivity */
    diffusivity (nx, ny, hx, hy, eps, u, dc);

    /* inner iterations: Gauss-Seidel */
    for (k=1; k<=maxgs; k++)
        gauss_seidel (alpha, nx, ny, hx, hy, g, dc, u);

    /* analyse image */
    analyse (u, nx, ny, &min, &max, &mean, &vari);
    printf("minimum:                 %8.2lf \n", min);
    printf("maximum:                 %8.2lf \n", max);
    printf("mean:                    %8.2lf \n", mean);
    printf("variance:                %8.2lf \n\n", vari);
    } /* for m */


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
/* open file and write header (incl. filter parameters) */
comment_line (comments, "# TV IMAGE RESTORATION WITH THE KACANOV METHOD\n");
comment_line (comments, "# outer iterations:  fixed point\n");
comment_line (comments, "# inner iterations:  Gauss-Seidel\n");
comment_line (comments, "# initial image:     %s\n", in);
comment_line (comments, "# alpha:             %8.2lf\n", alpha);
comment_line (comments, "# epsilon:           %8.6lf\n", eps);
comment_line (comments, "# FP iterations:     %8ld\n", maxfp);
comment_line (comments, "# GS iterations:     %8ld\n", maxgs);
comment_line (comments, "# minimum:           %8.2lf\n", min);
comment_line (comments, "# maximum:           %8.2lf\n", max);
comment_line (comments, "# mean:              %8.2lf\n", mean);
comment_line (comments, "# variance:          %8.2lf\n", vari);


/* write image */
write_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- disallocate storage ---- */

disalloc_matrix (g,  nx+2, ny+2);
disalloc_matrix (dc, nx+2, ny+2);
disalloc_matrix (u,  nx+2, ny+2);

return(0);
}

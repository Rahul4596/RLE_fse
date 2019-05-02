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
} /* dummies */

/*--------------------------------------------------------------------------*/

long dithering

  (float** mask,  /* mask (output) */
   float** u,     /* input image */
   long nx,       /* image x dimensions */
   long ny,       /* image y dimensions */
   long bx,       /* boundary size x */
   long by,       /* boundary size y */
   float hx,      /* spatial grid size x (for finite differences)*/
   float hy,      /* spatial grid size y (for finite differences)*/
   float density) /* desired pixel density in (0,1] */

/* Create dithered mask of Laplacian magnitude with Floyd-Steinberg. */
{
long i, j;          /* loop variables */
float** laplace;    /* array for Laplacian magnitude */
float max_laplace;  /* maximum Laplacien magnitude */
float avg_laplace;  /* average Laplacian magnitude */
float avg_rescaled; /* rescaled average value for density adaption */
float old, error;   /* auxiliary variables for  dithering */
float rx, ry;       /* auxiliary variables for Laplace computation */
long mask_points;   /* number of mask points after dithering */

alloc_matrix(&laplace, nx + 2 * bx, ny + 2 * by);

/* create dummy boundaries for u by mirroring */
dummies(u, nx, ny);

/* compute laplacian magnitude, maximum and average */
max_laplace = 0;
avg_laplace = 0;
rx = 1.0 / (hx * hx);
ry = 1.0 / (hy * hy);
for (j = by; j < ny + by; j++)
   for (i = bx; i < nx + bx; i++)
   {
      laplace[i][j] = fabs(
            rx * (u[i + 1][j] + u[i - 1][j]) + ry * (u[i][j + 1] + u[i][j - 1])
                  - 2.0 * (rx + ry) * u[i][j]);
      avg_laplace += laplace[i][j];
      if (laplace[i][j] > max_laplace)
      {
         max_laplace = laplace[i][j];
      }
   }
avg_laplace /= (float) (nx * ny);

printf("Computed squared Laplacian magnitude: avg %f, max %f\n",
      avg_laplace, max_laplace);

/* use a transformation of type x -> a*x with suitable a such that
 new average is density*255 */
avg_rescaled = 0.0;

for (i = bx; i < nx + bx; i++)
{
   for (j = by; j < ny + by; j++)
   {
      mask[i][j] = (density * 255.0) / avg_laplace * laplace[i][j];
      avg_rescaled += mask[i][j];
   }
}
avg_rescaled /= (float) (nx * ny);
printf("Average after rescaling: %f (%f*255=%f))\n",
      avg_rescaled, density, density * 255.0);

/* perform floyd-steinberg dithering */
mask_points = 0;
for (j = by; j < ny + by; j++)
{
   for (i = bx; i < nx + bx; i++)
   {

      old = mask[i][j];

      /* quantisation */
      if (mask[i][j] >= fabs(255.0 - mask[i][j]))
      {
         mask[i][j] = 255.0;
         mask_points++;
      }
      else
      {
         mask[i][j] = 0.0;
      }

      error = old - mask[i][j];

      /* error distribution */
      mask[i + 1][j] += 7.0 / 16.0 * error;
      mask[i][j + 1] += 5.0 / 16.0 * error;
      mask[i + 1][j + 1] += 1.0 / 16.0 * error;
      mask[i - 1][j + 1] += 3.0 / 16.0 * error;
   }
}

printf("created %ld mask points (desired: %ld)\n", mask_points,
      (long) roundf(density * nx * ny));

disalloc_matrix(laplace, nx + 2 * bx, ny + 2 * by);

return mask_points;
} /* dithering */

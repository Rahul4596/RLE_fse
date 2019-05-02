#include <getopt.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
long generate_mask(float** mask,
                  long nx, long ny,
                  long bx, long by,
                  long desired_pixels) {

  long i,j,x,y;

  /* initialise empty mask */
  for (j=by; j<ny+by; j++) {
    for (i=bx; i<nx+bx; i++) {
          mask[i][j] = 0.0;
    }
  }

  long current_pixels = 0;

  /* Add uniformly distributed random pixels using rand(). Make sure that
     you actually get the desired number of pixels. */
  while (current_pixels < desired_pixels) {
  
    x = rand() % nx + 1;
    y = rand() % ny + 1;
    
    if (mask[x][y] < 0.5) {
      mask[x][y] = 1.0;
      current_pixels++;
    }
  }

  return current_pixels;
}


/*--------------------------------------------------------------------------*/

/* Create regular mask that contains as many known data points as specified
in the variable "pixels" on a sampling grid. The grid distance in x- and y-
direction is adapted to the aspect ratio and aims to approximate the desired
pixel number as best as possible. */

long generate_regular_mask(float** mask,          /* mask (output) */
                           long nx, long ny,      /* image dimensions */
                           long bx, long by,      /* boundary size */
                           long pixels) {         /* desired number of pixels */

  long i,j;
  long hx,hy;
  long remainderx, remaindery;
  long pixelsx, pixelsy;

  /* initialise empty mask */
  for (j=by; j<ny+by; j++)
    for (i=bx; i<nx+bx; i++)
      {
        mask[i][j] = 0.0;
      }

  long current_pixels = 0;

  /* compute number of pixels per side: pixelsx * pixelsy should approximate
     pixels and the distance between known points should be similar in x- and
     y-direction  */
  pixelsx = (long)roundf(sqrtf((float)pixels*(float)nx/(float)ny));
  pixelsy = (long)roundf((float)pixelsx * (float)ny / (float) nx);
  printf("%ld pixels in x-direction, %ld pixels in y-direction (%ld/%ld)\n",
         pixelsx,pixelsy,pixelsx*pixelsy,pixels);
  
  /* compute distances hx/hy between pixels, e.g. hx=2 would give you 
     1 0 0 1 0 0 1 etc. Remember that the left/right/top/bottom points should
     be on the image boundary (i.e. coordinate bx/by and/or nx+bx-1/ny+by-1) */
  hx = (nx-pixelsx)/(pixelsx-1);
  hy = (ny-pixelsy)/(pixelsy-1);
  remainderx = (nx-pixelsx)%(pixelsx-1);
  remaindery = (ny-pixelsy)%(pixelsy-1);
  printf("%ld grid-size x (remainder %ld), %ld grid-size y (remainder %ld)\n",
         hx,remainderx,hy,remaindery);
  
  /* set mask points according to grid distances hx and hy. if 
     the exact distances cannot be replicated, use the remainders to enlarge
     as many grid distances as needed by 1. */
  for (i=bx; i<bx+nx; i+=hx+1) {
    for (j=by; j<by+ny; j+=hy+1)
      {
        mask[i][j] = 1.0;
        current_pixels++;
        if (remaindery>0) {
          remaindery--;
          j++;
        }
      }
    if (remainderx>0) {
      remainderx--;
      i++;
    }
    remaindery = (ny-pixelsy)%(pixelsy-1);
  }

  return current_pixels;
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
}  

/****************************************************************************/

void prepare_quantisation_map(long q, /* number of quantised values */
                              long max_q, /* maximum possible grey value */
                              long* quantisation_map_small,
                              long* quantisation_map_back) {
  long i;
  long a = 0;
  long max_levels=65536; /* 16 bit */

  //printf("preparing quantisation map with q=%ld and max_q=%ld\n",q,max_q); 
        
  if (q > max_levels) {
    printf("warning, quantisation exceeds max levels, set to %ld\n",max_levels);
    q = max_levels;
  }
        
  if ((q == max_levels || q == 0) && (max_q == max_levels)) {
      for (i = 0; i < q; ++i) {
        quantisation_map_back[i] = i;
        quantisation_map_small[i] = i;
      }
  }
  
  for (i = 0; i <= q; ++i) {
    quantisation_map_back[i] =
      (long)((i+0.5)*max_q/q);
  }
  //printf("map back: done\n");
        
  a=0;
  for (i = 0; i <= max_q-1; ++i) {
    if (i-quantisation_map_back[a] > quantisation_map_back[a+1]-i
        && a < q-1) a++;
    quantisation_map_small[i] = a;
  }
  //printf("map small: done\n");
}


/*--------------------------------------------------------------------------*/
void quantise_image_scalar(long** original,
                           long nx, long ny,
                           long max_q,
                           long q,
                           long* quantisation_map_small,
                           long* quantisation_map_back,
                           long** quantised) {
  
  long x, y;
  long bx = 1, by = 1;
  
  prepare_quantisation_map(q,max_q,quantisation_map_small,quantisation_map_back);

    for (x = bx; x < nx + bx; x++)
      for (y = by; y < ny + by; y++) {
        if (original[x][y]>255) {
          original[x][y]=255;
        }
        if (original[x][y]<0) {
          original[x][y]=0;
        }
        quantised[x][y] = quantisation_map_small[(long)(original[x][y])];
        /*if (quantised[c][x][y]>=q) {
          printf("%f > %ld\n",original[c][x][y],quantised[c][x][y]);
          }*/
      }
}


/*--------------------------------------------------------------------------*/
void init_image_quantised(long** quantised,
                          float** mask,
                          long nx, long ny,
                          long* qmap,
                          long** image) {
  
  long i, j;

  if (mask==0) {
    for (i=1;i<=nx;i++) {
      for (j=1;j<=ny;j++) {
          image[i][j]=(float)qmap[quantised[i][j]];
      }
    }
  } else {
    for (i=1;i<=nx;i++) {
      for (j=1;j<=ny;j++) {
        if (mask[i][j]>0.5) {
          image[i][j]=(float)qmap[quantised[i][j]];
        } else {
          image[i][j]=0;
        }
      }
    }
  }
}

/*--------------------------------------------------------------------------*/

/* Create dithered mask of Laplacian magnitude. */

long dithering(float** mask,          /* mask (output) */
               float** u,             /* input image */
               long nx, long ny,      /* image dimensions */
               long bx, long by,      /* boundary size */
               float hx, float hy,    /* spatial grid size (for finite differences)*/
               float density, long* qmap) {       /* desired pixel density in (0,1] */

  long i,j;
  float** laplace;
  float max_laplace;
  float avg_laplace;
  float avg_rescaled;
  float old,error;
  float rx,ry;
  long mask_points;
  
  alloc_matrix(&laplace, nx+2*bx, ny+2*by);

  /* create dummy boundaries for u by mirroring */
  dummies (u, nx, ny);
  
  /* compute laplacian magnitude, maximum and average */
  max_laplace = 0;
  avg_laplace = 0;
  rx = 1.0 / (hx * hx);
  ry = 1.0 / (hy * hy);
  for (j=by; j<ny+by; j++)
    for (i=bx; i<nx+bx; i++) {
      laplace[i][j]=fabs(rx*(u[i+1][j] + u[i-1][j]) + ry*(u[i][j+1] + u[i][j-1])
                         -2.0 *(rx+ry)*u[i][j]);
      avg_laplace+=laplace[i][j];
      if (laplace[i][j] > max_laplace) {
        max_laplace = laplace[i][j];
      }
    }
  avg_laplace/=(float)(nx*ny);

  printf("Computed squared Laplacian magnitude: avg %f, max %f\n",
         avg_laplace, max_laplace);

  /* use a transformation of type x -> a*x with suitable a such that
     new average is density*255 */
  avg_rescaled = 0.0;
  
  for (i=bx; i<nx+bx; i++) {
    for (j=by; j<ny+by; j++) {
      mask[i][j] = (density*255.0)/avg_laplace*laplace[i][j];
      avg_rescaled+=mask[i][j];
      // printf("old: %f, rescaled: %f\n",laplace[i][j],mask[i][j]);
    }
  }
  avg_rescaled/=(float)(nx*ny);
  printf("Average after rescaling: %f (%f*255=%f))\n",
     avg_rescaled,density,density*255.0);

  /* perform floyd-steinberg dithering */
  mask_points = 0;
  long ref[2], run, flag, k;
  ref[0]=bx;
  ref[1]=by;
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
         /*flag=0;
         run=(j-ref[1])*nx + (i-ref[0]);
         for(k=0;k<40;k++)
         {
          if(run==qmap[k])
            {
              
              ref[0]=i;
              ref[1]=j;
              break;
            }
         }*/

         
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
         (long)roundf(density*nx*ny));

  disalloc_matrix(laplace, nx+2*bx, ny+2*by);

  return mask_points;
}

/*--------------------------------------------------------------------------*/

void comment_line

     (char* comment,       /* comment string (output) */
      char* lineformat,    /* format string for comment line */
      ...)                 /* optional arguments */

/* 
  Add a line to the comment string comment. The string line can contain plain
  text and format characters that are compatible with sprintf.
  Example call: print_comment_line(comment,"Text %f %ld",float_var,int_var);
  If no line break is supplied at the end of the input string, it is added
  automatically.
*/

{
char     line[80];
va_list  arguments;

/* get list of optional function arguments */
va_start (arguments, lineformat);

/* convert format string and arguments to plain text line string */
vsprintf (line, lineformat, arguments);

/* add line to total commentary string */
strncat (comment, line, 80);

/* add line break if input string does not end with one */
if (line[strlen(line)-1] != '\n')
   sprintf (comment, "%s\n", comment);

/* close argument list */
va_end (arguments);

return;

} /* comment_line */


/*--------------------------------------------------------------------------*/

void write_mask

(float  **u,          /* image, unchanged */
 long   nx,           /* image size in x direction */
 long   ny,           /* image size in y direction */
 long   known,        /* which value should be assigned to known pixels? */
 char   *file_name,   /* name of pgm file */
 char   *comments)    /* comment string (set 0 for no comments) */

/*
  writes a binary mask image into a pgm P5 file;
 */

{
  FILE           *outimage;  /* output file */
  long           i, j;       /* loop variables */
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
        if (known > 0) {
          if (u[i][j] < 0.5)
            byte = (unsigned char)(255.0);
          else
            byte = (unsigned char)(0.0);
        } else {
          if (u[i][j] < 0.5)
            byte = (unsigned char)(0.0);
          else
            byte = (unsigned char)(255.0);
        }
        fwrite (&byte, sizeof(unsigned char), 1, outimage);
      }

  /* close file */
  fclose (outimage);
  
  return;
  
} /* write_pgm */

/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

int main(int argc, char** args) {
  int ch;	  /* all-purpose char */
  char used[256]; /* flag for all possible chars, used to check for multiple
                   * occurences of input parameters */
  long nx=256,ny=256; /* image dimensions */
  float density = 0.1;
  long desired_pixels;
  float** mask=0;
  float** image=0;
  char* image_file=0, *mask_file = 0;
  long bx =1, by =1;
  char comments[256];
  long i,j;
  long regular_flag = 0;
  long invert_flag = 0;
  long floyd_steinberg = 0;
  
  /* process input parameters */
  for (ch = 0;ch <= 255; ch++)
    used[ch] = 0;
  while ((ch = getopt(argc,args,"ri:m:x:y:d:If")) != -1) {
    used[ch]++;

    if (used[ch] > 1) {
      printf("Double parameter: %c\n",ch);
      printf("Please check your input again.\n");
    }

    switch(ch) {
      case 'i': image_file = optarg;break;
      case 'm': mask_file = optarg;break;
      case 'r': regular_flag = 1;break;
      case 'x': nx = atoi(optarg);break;
      case 'y': ny = atoi(optarg);break;
      case 'd': density = atof(optarg);break;
      case 'f': floyd_steinberg = 1;break;
      case 'I': invert_flag = 1; break;
      default: printf("Unknown argument.\n");
    }
  }

  if (mask_file==0) {
    printf("Mask file must be specified\n");
    return 0;
  }

  /* load images */
  if (image_file != 0) {
    read_pgm_and_allocate_memory (image_file, &nx, &ny, &image);
  }

  desired_pixels = (long)round(density*nx*ny);

  alloc_matrix(&mask, nx+2*bx, ny+2*by);

  for (j=by; j<ny+by; j++)
    for (i=bx; i<nx+bx; i++) {
      mask[i][j]=0.0;
    }

    long qmap[257], qmap_small[40];
  prepare_quantisation_map(40,256,qmap,qmap_small);

  if (floyd_steinberg == 1) {
    dithering(mask,image,nx,ny,bx,by,1.0,1.0,density,qmap_small);
    if (invert_flag>0) invert_flag = 0; else invert_flag = 1; 
  } else if (regular_flag == 0) {
    generate_mask(mask,nx,ny,bx,by,desired_pixels);
  } else {
    printf("regular mode\n");
    generate_regular_mask(mask,nx,ny,bx,by,desired_pixels);
  }

  //long qmap[257], qmap_small[20];
  //prepare_quantisation_map(20,256,qmap,qmap_small);

  //int i;
  for(i=0;i<20;i++)
  	printf("%ld -> %ld \n",i,qmap_small[i] );

  comments[0]='\0';
  comment_line(comments,"#mask with density %f, regular flag %ld",
               density,regular_flag);
  write_mask(mask,nx,ny,invert_flag,mask_file,comments);
  printf ("output mask %s successfully written\n\n", mask_file);

  

  
  return 0;
}

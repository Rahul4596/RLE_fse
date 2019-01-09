/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                  PROGRAMMING EXERCISE FOR THE LECTURE                    */
/*                            IMAGE COMPRESSION                             */
/*                               JPEG LIGHT                                 */
/*                  (Copyright by Pascal Peter, 6/2017)                     */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* global includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <getopt.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

/* local includes */
#include "alloc.h"              /* memory allocation */
#include "image_io.h"           /* reading and writing pgm and ppm images */
#include "bfio.h"               /* writing and reading of bitfiles */

/* defines */
/* version */
#define VERSION 1.1-06-16
/* supported input formats */
#define FORMAT_PGM 0
#define FORMAT_PPM 1
/* auxiliary functions */
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
/* maximum number of channels */
#define MAXCHANNELS 3
/* define maximum grey value */
#define MAXGREYVALUE 255
/* console formatting */
#define ONE_UP "\033[5D\033[1A"
#define CLRLINE "\033[K"
/* symbols */
#define ZRL 192
#define EOB 193

/* definition of compressed image datatype and struct */
typedef struct ImageData ImageData;
struct ImageData {
  long*** orig_rgb;       /* original image (RGB) */
  long*** orig_ycbcr;     /* original image (YCbCr) */
  long size_orig;         /* size of raw input image */
  long nx, ny, nc;        /* image dimensions and channels */
  long nx_ext[3], ny_ext[3]; /* extended image sizes */
  long block_size;        /* size of quadratic DCT blocks */
  long blocks_x,          /* number of blocks in each direction */
       blocks_y;
  double*** dct;          /* DCT coefficients of original image */
  long***   dct_quant;    /* quantised DCT coefficients */
  double*** rec;          /* reconstruction from DCT  */
  long***   rec_quant;    /* integer reconstruction  */
  long      s;            /* subsampling factor */
};

long log2long(long x) {
  int logx = 0;
  while (x >>= 1) ++logx;
  return logx;
}

/*--------------------------------------------------------------------------*/
void init_image (ImageData* img) {
  /* set image variables to default values */
  img->nx=img->ny=img->nc=0;
  /* img->nx_ext */
  /*   =img->ny_ext=0; */
  img->blocks_x=0;img->blocks_y=0;
  img->size_orig=0;
  img->orig_rgb=0;
  img->orig_ycbcr=0;
  img->rec=0;
  img->rec_quant=0;
  img->dct=0;
  img->dct_quant=0;
  img->block_size=8;
  img->s=2;
}

/*--------------------------------------------------------------------------*/
void alloc_image (ImageData* img,long nx,long ny) {

  long blocks_x, blocks_y, nx_ext, ny_ext, N, nx_coarse, ny_coarse;

  /* set number of blocks and extended image size */
  N = img->block_size;
  nx = img->nx;
  ny = img->ny;
  blocks_x = nx/img->block_size;
  if ((nx % N) > 0) blocks_x++;
  blocks_y = ny/img->block_size;
  if ((ny % N) > 0) blocks_y++;
  nx_ext = blocks_x*img->block_size;
  ny_ext = blocks_y*img->block_size;

  img->nx_ext[0] = nx_ext;
  img->ny_ext[0] = ny_ext;
  img->blocks_x = blocks_x;
  img->blocks_y = blocks_y;

  
  /* determine coarse resolution */
  nx_coarse = nx/img->s;
  ny_coarse = ny/img->s;
  if ((nx % img->s) > 0) nx_coarse++;
  if ((ny % img->s) > 0) ny_coarse++;

  /* determine extended coarse resolution */
  blocks_x = nx_coarse/img->block_size;
  if ((nx_coarse % N) > 0) blocks_x++;
  blocks_y = ny_coarse/img->block_size;
  if ((ny_coarse % N) > 0) blocks_y++;
  img->nx_ext[1] = blocks_x*img->block_size;
  img->ny_ext[1] = blocks_y*img->block_size;
  img->nx_ext[2] = img->nx_ext[1];
  img->ny_ext[2] = img->ny_ext[1];

  /* printf("nx_ext %ld %ld %ld, ny_ext %ld %ld %ld\n", */
  /*        img->nx_ext[0],img->nx_ext[1],img->nx_ext[2],img->ny_ext[0], */
  /*        img->ny_ext[1],img->ny_ext[2]); */
  
  /* allocate all image data arrays */
  if (img->orig_rgb == 0)
    alloc_long_cubix(&img->orig_rgb,MAXCHANNELS,nx_ext+2,ny_ext+2);
  if (img->orig_ycbcr == 0)
    alloc_long_cubix(&img->orig_ycbcr,MAXCHANNELS,nx_ext+2,ny_ext+2);
  if (img->dct_quant == 0)
    alloc_long_cubix(&img->dct_quant,MAXCHANNELS,nx_ext+2,ny_ext+2);
  if (img->dct == 0)
    alloc_cubix(&img->dct,MAXCHANNELS,nx_ext+2,ny_ext+2);
  if (img->rec == 0)
    alloc_cubix(&img->rec,MAXCHANNELS,nx_ext+2,ny_ext+2);
  if (img->rec_quant == 0)
    alloc_long_cubix(&img->rec_quant,MAXCHANNELS,nx_ext+2,ny_ext+2);
}

/*--------------------------------------------------------------------------*/
void destroy_image (ImageData* img) {
  /* disalloc all image data arrays */
  long nx_ext, ny_ext;
  nx_ext = img->nx;
  ny_ext = img->ny;

  /* disalloc images */
  if (img->orig_rgb != 0)
    disalloc_long_cubix(img->orig_rgb,MAXCHANNELS,nx_ext+2,ny_ext+2);
  if (img->orig_ycbcr != 0)
    disalloc_long_cubix(img->orig_ycbcr,MAXCHANNELS,nx_ext+2,ny_ext+2);
  if (img->dct != 0)
    disalloc_cubix(img->dct,MAXCHANNELS,nx_ext+2,ny_ext+2);
  if (img->dct_quant != 0)
    disalloc_long_cubix(img->dct_quant,MAXCHANNELS,nx_ext+2,ny_ext+2);
  if (img->rec != 0)
    disalloc_cubix(img->rec,MAXCHANNELS,nx_ext+2,ny_ext+2);
  if (img->rec_quant != 0)
    disalloc_long_cubix(img->rec_quant,MAXCHANNELS,nx_ext+2,ny_ext+2);
}

/*--------------------------------------------------------------------------*/
void write_comment_string(ImageData* img, char* additional, char* comments)
/* writes all important information of an R-EED image into a comment string
 * parameter additional allows to add more custom text, set to 0 if not needed.
 * make sure that the comments array is preallocated with sufficient space. */
{
  /* write parameter values in comment string */
  comments[0] = '\0';
  comment_line (comments, (char*)"# IMAGE COMPRESSION - PROGRAMMING EXERCISE\n");
}

/*--------------------------------------------------------------------------*/
/* for user interaction */
/* Returns true if it is possible to read from cin */
int canread()
{
  struct timeval tv;
  tv.tv_sec = 0;
  tv.tv_usec = 0;
  fd_set read_set;
  FD_ZERO(&read_set);
  FD_SET(0,&read_set);
  select(1,&read_set,NULL,NULL, &tv);
  return (FD_ISSET(0,&read_set));
}

/*--------------------------------------------------------------------------*/
void copy_cubix(double ***source, double ***target,
    long sc, long nc, long nx, long ny) {
  /* copy cubix from channel sc to channel nc */
  long i,j,m;
  for (m=0;m<nc;m++) {
    for (i=1;i<=nx;i++) {
      for (j=1;j<=ny;j++) {
        target[m][i][j]=source[m][i][j];
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void copy_cubix_long(long ***source, long ***target,
                    long sc, long nc, long nx, long ny) {
  /* copy cubix (only nc channels, starting with channel sc) */
  long i,j,m;
  for (m=sc;m<nc;m++) {
    for (i=1;i<=nx;i++) {
      for (j=1;j<=ny;j++) {
        target[m][i][j]=source[m][i][j];
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void copy_matrix_long(long **source, long **target,long nx, long ny) {
  /* copy input matrix to target matrix */
  long i,j;
  for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
      target[i][j]=source[i][j];
    }
  }
}

/*--------------------------------------------------------------------------*/
void copy_vector_long(long *source, long *target, long nx) {
  /* copy input vector to target vector */
  long i;
  for (i=0;i<nx;i++) {
    target[i]=source[i];
  }
}

/*--------------------------------------------------------------------------*/
void convert_matrix_int(double **source, long **target,long nx, long ny) {
  /* copy input matrix to target matrix */
  long i,j;
  for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
      target[i][j]=(long)round(source[i][j]);
    }
  }
}

/*--------------------------------------------------------------------------*/
long get_size_of_file(char* filename) {
  /* compute and return size of file with arbitrary content */
  FILE   *infile;             /* input file */
  long i;
  /* open file */
  infile = fopen(filename,"r");
  if (NULL == infile)
    {
      printf("Could not open file '%s' for reading, aborting (10)\n",filename);
      exit(-1);
    }
  /* Find length of file */
  fseek(infile,0,SEEK_END);
  i = ftell(infile);
  fclose(infile);
  return i;
}

/*--------------------------------------------------------------------------*/
double get_compression_ratio(char* original_file,
                             char* compressed_file) {
  /* computes compression ratio between two files on the hard disk (X:1) */
  long original_size;
  long compressed_size;
  
  original_size = get_size_of_file(original_file);
  compressed_size = get_size_of_file(compressed_file);
  
  return (double)original_size/(double)compressed_size;
}

/*--------------------------------------------------------------------------*/
void convert_image_to_vector (long*** image_matrix,
                              long bx, long by, long nx, long ny, long nc,
                              long* image_vector) {
  long i,j,c;

  /* parse image channel by channel, row by row */
  for (c = 0; c < nc; c++) 
    for (j = by; j < ny+by;j++)
      for (i = bx; i < nx+bx;i++) {
        image_vector[c*nx*ny+(j-by)*nx+(i-bx)]=image_matrix[c][i][j];
      }
}

/*--------------------------------------------------------------------------*/
void convert_vector_to_image (long* image_vector,
                              long bx, long by, long nx, long ny, long nc,
                              long*** image_matrix) {
  long i,j,c;

  /* parse image channel by channel, row by row */
  for (c = 0; c < nc; c++) 
    for (j = by; j < ny+by;j++)
      for (i = bx; i < nx+bx;i++) {
        image_matrix[c][i][j]=image_vector[c*nx*ny+(j-by)*nx+(i-bx)];
        /*  printf("%ld %ld %ld (%ld): %ld\n",
           c,i,j,c*nx*ny+(j-by)*nx+(i-bx),image_matrix[c][i][j]); */ 
      }
}

/*--------------------------------------------------------------------------*/
void abs_img(long** image,long nx, long ny, long** normalised) {
  long i,j;

  for (i=1;i<=nx;i++)
    for (j=1;j<=ny;j++) {
      normalised[i][j]=abs(image[i][j]);
    }
}

/*--------------------------------------------------------------------------*/
void normalise_to_8bit(long** image,long nx, long ny, long** normalised) {
  long i,j;
  long min, max;

  min=1000000;
  max=-1000000;
  for (i=1;i<=nx;i++)
    for (j=1;j<=ny;j++) {
      if (min>image[i][j]) min = image[i][j];
      if (max<image[i][j]) max = image[i][j];
    }

  for (i=1;i<=nx;i++)
    for (j=1;j<=ny;j++) {
      normalised[i][j]=((image[i][j]-min)*255)/(max-min);
    }
}

/*--------------------------------------------------------------------------*/
void print_usage_message() {
  printf("./compress -i input_file -o output_prefix [optional parameters]\n");
  printf("list of mandatory paramters:\n");
  printf("-i input_file   (string): uncompressed image, e.g. \"image.ppm\"\n");
  printf("-o out_prefix   (string): prefix for output files, e.g. \"my_image\"\n");
  printf("list of optional paramters:\n");
  printf("-s subsampling factor      (int): subsample each dimension by factor s\n");
  printf("-q quantisation parameter  (int): use uniform quantisation matrix with entry q everywhere\n");
}

/*--------------------------------------------------------------------------*/

/* apply WNC algorithm for adaptive arithmetic integer encoding */
void encode_adaptive_wnc(
  long* sourceword,   /* array containing n numbers from {0,...,s}
                         where s is the end of file symbol */
  long n,             /* length of sourceword */
  long s,             /* size of source alphabet */
  float r,           /* rescaling parameter */
  long  M,            /* WNC discretisation parameter M */
  FILE* debug_file,   /* 0 - no output, otherwise debug output to file */
  BFILE* compressed)  {/* binary file for compressed bitstring */

  long i,j;      /* loop variables */
  long L,H;      /* low and high interval endpoints of current interval */
  long oldL;     /* temporary variable to preserve L for computing new interval*/
  long C;        /* sum of all counters */
  long k;        /* underflow counter */
  long* counter; /* array of counters for adaptive probabilities */
  long M12=M/2, M14=M/4, M34=3*M/4;  /* time savers */
  long symbol;   /* index of current symbol in counter array */
  long csum;     /* sum of counters 0,...,symbol-1 */
  long rescaling_counter = 0;
  long bits_written = 0;


  /* allocate memory */
  alloc_long_vector(&counter,s);

  /* initialise counters and C */
  for (i=0;i<s;i++) {
    counter[i]=1;
  }
  C = s;

  if (debug_file != 0) {
    fprintf(debug_file,"n: %ld, s: %ld, r: %f, M: %ld\n",n,s,r,M);
  }

  if ((double)C>(double)M/4.0+2.0) {
    if (debug_file != 0) {
      fprintf(debug_file,
              "M=%ld is too small (C=%ld), setting M to %ld\n",M,C,8*C);
    }
    M=8*C;
    M12=M/2; M14=M/4; M34=3*M/4;
  }

  /* initialise interval endpoints and underflow counter */
  L=0;
  H=M;
  k=0;

  /* encode sourceword */
  for (i=0;i<n;i++) {
    if (debug_file != 0) {
      fprintf(debug_file,"sourceword[%ld]=%ld\n",i,sourceword[i]);
    }
    /* underflow expansions/rescaling */
    while (1) {
      /* check for underflow expansion x -> 2*x - M/2 */
      if ((L >= M14) && (L<M12) && (H>M12) && (H<=M34)) {
        L=2*L-M12; H=2*H-M12;
        k++;
        if (debug_file != 0) {
          fprintf(debug_file,"underflow: x -> 2*x - %ld\n",M12);
          rescaling_counter++;
        }
        continue;
      }

      /* check for rescaling x -> 2*x */
      if (H<=M12) {
        if (debug_file != 0) {
          fprintf(debug_file,"rescaling: x-> 2*x:");
          rescaling_counter++;
        }
        L=2*L; H=2*H;
        /* write 01^k to bitstream */
        set_bit(compressed,0);
        if (debug_file != 0) {
          fprintf(debug_file,"written bits: 0");
        }
        for (j=0;j<k;j++) {
          set_bit(compressed,1);
          if (debug_file != 0) {
            fprintf(debug_file,"1");
          }
        }
        if (debug_file != 0) {
          fprintf(debug_file,"\n");
        }
        k=0;
        continue;
      }

      /* check for rescaling x -> 2*x - M/2 */
      if (L>=M12) {
        if (debug_file != 0) {
          fprintf(debug_file,"rescaling: x-> 2*x - %ld:",M);
          rescaling_counter++;
        }
        L=2*L-M; H=2*H-M;
        /* write 10^k to bitstream */
        set_bit(compressed,1);
        if (debug_file != 0) {
          fprintf(debug_file,"written bits: 1");
        }
        for (j=0;j<k;j++) {
          set_bit(compressed,0);
          if (debug_file != 0) {
            fprintf(debug_file,"0");
          }
        }
        if (debug_file != 0) {
          fprintf(debug_file,"\n");
        }
        k=0;
        continue;
      }

      if (debug_file != 0) {
        fprintf(debug_file,"k: %ld, [%ld, %ld)\n",k,L,H);
      }
      break;
    }

    /* readjustment */
    while ((double)C>(double)M/4.0+2.0) {
      if (debug_file != 0) {
        fprintf(debug_file,"C: %ld, M/4+2.0=%f\n",C,M/4.0+2.0);
      }
      C=0;
      for (j=0;j<s;j++) {
        counter[j]=(long)round((double)counter[j]*r);
        if (counter[j]==0) counter[j]=1; /* no zero-counters allowed */
        C+=counter[j];
      }
    }

    /* encode symbol */
    symbol=sourceword[i];
    csum=0;
    for (j=0;j<symbol;j++) csum+=counter[j];

    oldL=L;
    L=L+(long)floor((double)(csum*(H-L))/(double)C);

H=oldL+(long)floor((double)((csum+counter[symbol])*(H-oldL))/(double)C);
    if (debug_file != 0) {
      fprintf(debug_file,"new [L,H) = [%ld,%ld)\n",L,H);
    }
    counter[symbol]++; C++;
  }

      /* underflow expansions/rescaling */
    while (1) {
      /* check for underflow expansion x -> 2*x - M/2 */
      if ((L >= M14) && (L<M12) && (H>M12) && (H<=M34)) {
        L=2*L-M12; H=2*H-M12;
        k++;
        if (debug_file != 0) {
          fprintf(debug_file,"underflow: x -> 2*x - %ld\n",M12);
        }
        continue;
      }

      /* check for rescaling x -> 2*x */
      if (H<=M12) {
        if (debug_file != 0) {
          fprintf(debug_file,"rescaling: x-> 2*x:");
        }
        L=2*L; H=2*H;
        /* write 01^k to bitstream */
        set_bit(compressed,0);
        bits_written++;
        if (debug_file != 0) {
          fprintf(debug_file,"written bits: 0");
        }
        for (j=0;j<k;j++) {
          set_bit(compressed,1);
          bits_written++;
          if (debug_file != 0) {
            fprintf(debug_file,"1");
          }
        }
        if (debug_file != 0) {
          fprintf(debug_file,"\n");
        }
        k=0;
        continue;
      }

      /* check for rescaling x -> 2*x - M/2 */
      if (L>=M12) {
        if (debug_file != 0) {
          fprintf(debug_file,"rescaling: x-> 2*x - %ld:",M);
        }
        L=2*L-M; H=2*H-M;
        /* write 10^k to bitstream */
        set_bit(compressed,1);
        bits_written++;
        if (debug_file != 0) {
          fprintf(debug_file,"written bits: 1");
        }
        for (j=0;j<k;j++) {
          set_bit(compressed,0);
          bits_written++;
          if (debug_file != 0) {
            fprintf(debug_file,"0");
          }
        }
        if (debug_file != 0) {
          fprintf(debug_file,"\n");
        }
        k=0;
        continue;
      }

      if (debug_file != 0) {
        fprintf(debug_file,"k: %ld, [%ld, %ld)\n",k,L,H);
      }
      break;
    }


  /* last step */
  if (debug_file != 0) {
    fprintf(debug_file,"last interval - written bits:");
  }
  if (L<M14) {
      set_bit(compressed,0);
      bits_written++;
      if (debug_file != 0) {
        fprintf(debug_file,"0");
      }
    for (j=0;j<k+1;j++) {
      set_bit(compressed,1);
      bits_written++;
      if (debug_file != 0) {
        fprintf(debug_file,"1");
      }
    }
  } else {
    set_bit(compressed,1);
    bits_written++;
    if (debug_file != 0) {
      fprintf(debug_file,"1");
    }
    for (j=0;j<k+1;j++) {
      set_bit(compressed,0);
      bits_written++;
      if (debug_file != 0) {
        fprintf(debug_file,"0");
      }
    }
  }

  /* write additional bits that have to be read by decoder */
  if (debug_file != 0) {
    fprintf(debug_file,"\n padding bits:");
  }
  for (i=0;i<log2(M)-bits_written;i++) {
    set_bit(compressed,0);
    if (debug_file != 0) {
      fprintf(debug_file,"0");
    }
  }

  if (debug_file != 0) {
    fprintf(debug_file,"\n rescalings/underflow expanions: %ld, additional bits_written: %ld\n",
            rescaling_counter,bits_written);
  }

  /* allocate memory */
  disalloc_long_vector(counter,s);
}



/*--------------------------------------------------------------------------*/

/* apply WNC algorithm for adaptive arithmetic integer encoding */
void decode_adaptive_wnc(
    BFILE* compressed,  /* binary file with compressed bitstring */
    long n,             /* length of sourceword */
    long s,             /* size of source alphabet */
    double r,           /* rescaling parameter */
    long  M,            /* WNC discretisation parameter M */
    FILE* debug_file,   /* 0 - no output, 1 - debug output to file */
    long* sourceword)  {/* array containing n numbers from {0,...,s}
                           where s is the end of file symbol */
  long i,j;      /* loop variables */
  long L,H;      /* low and high interval endpoints of current interval */
  long oldL;     /* temporary variable to preserve L for computing new interval*/
  long C;        /* sum of all counters */
  long* counter; /* array of counters for adaptive probabilities */
  long M12=M/2, M14=M/4, M34=3*M/4;  /* time savers */
  long symbol;   /* index of current symbol in counter array */
  long csum;     /* sum of counters 0,...,symbol-1 */
  long w;        /* variable for finding correct decoding inverval */
  long v;        /* partial dyadic fraction */
  long b;        /* auxiliary variable for reading individual bits */
  long N;        /* auxiliary variable for determining number of initial bits
                    for v */
 
  printf("n=%ld s=%ld M=%ld \n",n,s,M );
  /* allocate memory */
  alloc_long_vector(&counter,s);
  
  /* initialise counters and C */
  for (i=0;i<s;i++) {
    counter[i]=1;
  }
  C = s;

  if ((double)C>(double)M/4.0+2.0) {
    if (debug_file != 0) {
      fprintf(debug_file,
              "M=%ld is too small (C=%ld), setting M to %ld\n",M,C,8*C);
    }
    M=8*C;
    M12=M/2; M14=M/4; M34=3*M/4;
  }

  /* initialise interval endpoints*/
  L=0;
  H=M;

  /* read first bits of codeword to obtain initival v */
  N=log2long(M); /* assumes that M is a power of 2! */
  j=(long)pow(2,N-1);
  v=0;
  for (i=0;i<N;i++) {
   b = get_bit(compressed);
   v += b*j;
   if (debug_file != 0) {
     fprintf(debug_file,
             "v: %ld, j: %ld, i: %ld, b: %ld\n",v,j,i,b);
   }
   j/=2;
  }
  if (debug_file != 0) {
    fprintf(debug_file,"initial v: %ld (%ld first bits from coded file, %f)\n",
           v,N,log((double)M)/log(2.0));
  }

  /* decode sourceword */
  for (i=0;i<n;i++) {
  
    /* underflow expansions/rescaling */
    while (1) {
      /* check for underflow expansion x -> 2*x - M/2 */
      if ((L >= M14) && (L<M12) && (H>M12) && (H<=M34)) {
        L=2*L-M12; H=2*H-M12; v=2*v-M12;
        /* shift in next bit */
        b=get_bit(compressed);
        if (b!=EOF) {
          v+=b;
        }
        if (debug_file != 0) {
          fprintf(debug_file,
                 "underflow: x -> 2*x - %ld, [%ld,%ld), b %ld v %ld\n",
                 M12,L,H,b,v);
        }
        continue;
      }
    
      /* check for rescaling x -> 2*x */
      if (H<=M12) {
        L=2*L; H=2*H; v=2*v;
        /* shift in next bit */
        b=get_bit(compressed);
        if (b!=EOF) {
          v+=b;
        }
        if (debug_file != 0) {
          fprintf(debug_file,"rescaling: x-> 2*x, [%ld, %ld), b %ld, v %ld\n",
                  L,H,b,v);
        }
        continue;
      }
    
      /* check for rescaling x -> 2*x - M/2 */
      if (L>=M12) {
        L=2*L-M; H=2*H-M; v=2*v-M;
        /* shift in next bit */
        b=get_bit(compressed);
        if (b!=EOF) {
          v+=b;
        }
        if (debug_file != 0) {
          fprintf(debug_file,
                  "rescaling: x-> 2*x - %ld, [%ld,%ld), b %ld, v %ld\n",
                 M,L,H,b,v);
        }
        continue;
      }

      if (debug_file != 0) {
        fprintf(debug_file,"v: %ld, [%ld, %ld)\n",v,L,H);
      }
      break;
    }
    
    /* readjustment */
    while ((double)C>(double)M/4.0+2.0) {
      if (debug_file != 0) {
        fprintf(debug_file,"readjust C: %ld, M/4+2.0=%f\n",C,M/4.0+2.0);
      }
      C=0;
      for (j=0;j<s;j++) {
        counter[j]=(long)round((double)counter[j]*r);
        if (counter[j]==0) counter[j]=1; /* no zero-counters allowed */
        C+=counter[j];
      }
    }

    /* decode symbol */
    w=((v-L+1)*C-1)/(H-L);

    /* find correct interval */
    symbol=0;
    csum=0;
    while ((symbol < s-1) && ((csum > w) || (csum+counter[symbol])<=w)) {
      csum+=counter[symbol];
      symbol++;
    }
    sourceword[i]=symbol;
    oldL=L;
    L=L+(long)floor((double)(csum*(H-L))/(double)C);
    H=oldL+(long)floor((double)((csum+counter[symbol])*(H-oldL))/(double)C);
    counter[symbol]++; C++;
    if (debug_file != 0) {
      fprintf(debug_file,"[c_i,c_i-1) = [%ld %ld) ",csum,csum+counter[symbol]);
      fprintf(debug_file,"w: %ld symbol[%ld]: %ld, new [L,H)=[%ld,%ld)\n",
              w,i,symbol,L,H);
    }
  }
}


/*--------------------------------------------------------------------------*/
void RGB_to_YCbCr(long ***rgb, long ***ycbcr,long nx, long ny) {
  /* convert with modified YUV conversion formula of JPEG2000 */
  long tmp0, tmp1; /* temporary variable that avoids problems if in and out
                       array are identical */
  long i,j;
  for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
      tmp0 = rgb[0][i][j];
      tmp1 = rgb[1][i][j];
      ycbcr[0][i][j]=(rgb[0][i][j]+2*rgb[1][i][j]+rgb[2][i][j])/4;
      ycbcr[1][i][j]= rgb[2][i][j]-tmp1;
      ycbcr[2][i][j]= tmp0-tmp1;
    }
  }
}

/*--------------------------------------------------------------------------*/
void YCbCr_to_RGB(long ***ycbcr, long ***rgb,long nx, long ny) {
  /* convert with modified YUV conversion formula of JPEG2000 */
  long tmp; /* temporary variable that avoids problems if in and out array
               are identical */
  long i,j;
  for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
      tmp = ycbcr[1][i][j];
      rgb[1][i][j]=ycbcr[0][i][j]-(ycbcr[1][i][j]+ycbcr[2][i][j])/4;
      rgb[0][i][j]=ycbcr[2][i][j]+rgb[1][i][j];
      rgb[2][i][j]=tmp+rgb[1][i][j];
    }
  }
}

/*--------------------------------------------------------------------------*/
/* quantisation matrix */

long w[8][8]= {{10,15,25,37,51,66,82,100}, 
               {15,19,28,39,52,67,83,101},
               {25,28,35,45,58,72,88,105},
               {37,39,45,54,66,79,94,111},
               {51,52,58,66,76,89,103,119},
               {66,67,72,79,89,101,114,130},
               {82,83,88,94,104,114,127,142},
               {100,101,105,111,119,130,142,156}};


/*--------------------------------------------------------------------------*/

void extend_image(long **orig,               /* original input image */
                  long nx, long ny,          /* original image dimensions */
                  long N,                    /* size of quadratic blocks */
                  long **img_ext)            /* extended image */
{
  long x,y,k;              /* loop variables */
  long blocks_x, blocks_y; /* number of blocks in each direction */
  
  /* determine number of blocks */
  blocks_x = nx/N;
  if ((nx % N) > 0) blocks_x++;
  blocks_y = ny/N;
  if ((ny % N) > 0) blocks_y++;
  
  /* copy image into extended array */
  for (x=1; x<=nx; x++) {
    for (y=1; y<=ny; y++) {
      img_ext[x][y]=orig[x][y];
    }
  }

  /* extend temporary array in x-direction if necessary by mirroring,
     e.g. for a,b,c known, continue image by c,b,a */
  for (x=nx+1;x<=blocks_x*N;x++) {
    k=0;
    for (y=1;y<=ny;y++) {
      img_ext[x][y]=orig[nx-k][y];
    }
    k++;
  }

  /* extend temporary array in y-direction if necessary by mirroring,
     e.g. for a,b,c known, continue image by c,b,a */
  for (y=ny+1;y<=blocks_y*N;y++) {
    k=0;
    for (x=1;x<=blocks_x*N;x++) {
      img_ext[x][y]=orig[x][ny-k];
    }
    k++;
  }
}







/*--------------------------------------------------------------------------*/

/* write number in n-bit notation */
void write_long_bitwise(long c, /* (positive) number to write */
                        long n, /* number of bits */
                        FILE* debug_file,   /* 0 - no output, 
                                            otherwise debug output to file */
                        BFILE* output_file) /* 0 no output,
                                               otherwise write to binary file */
{
  long i;
  if (debug_file !=0) {
    for (i=0;i<n;i++) {
      fprintf(debug_file,"%ld",c%2);
      c/=2;
    }
  }
  fprintf(debug_file," ");
  if (output_file !=0) {
    for (i=0;i<n;i++) {
      set_bit(output_file,c%2);
      c/=2;
    }
  }
  return;
}

/*--------------------------------------------------------------------------*/
void read_long_bitwise(long *c, /* (positive) number to write */
                        long n, /* number of bits */
                        FILE* debug_file,   /* 0 - no output, 
                                            otherwise debug output to file */
                        BFILE* input_file)
{
  long temp=0,i=0;
  long t;
  for(i=0;i<n;i++)
  {
    t=bfgetb(input_file);
    fprintf(debug_file, "%ld",t );
    //temp+=(long)(pow(2.0,i))*t;
  }
  *c=temp;
  //if(debug_file!=0)
    
}
/*--------------------------------------------------------------------------*/
/* calulates block DCT of input image/channel */
void RLE_encode(long  **image,      /* input quantised DCT coefficients */
                  long nx, long ny,   /* image dimensions */
                  FILE* debug_file,
                  BFILE* fp)  /* 0 - no output, 
                                         otherwise debug output to file */
                //  BFILE *binary_file) 
                {/* file for binary output */

  long u,v,k,l,j,i=0;      /* loop variables */
  long N=8;            /* block size */
  long blocks_x;       /* number of blocks in each direction */
  long blocks_y;
  long ox,oy;          /* block offsets */
  long category_lookup[2048]; /* lookup table for categories */
  long zigzag_x[64];   /* x-index for zig-zag traversal of blocks */ 
  long zigzag_y[64];   /* y-index for zig-zag traversal of blocks */
  long last_dc = 0;    /* previously encoded dc coefficient */
  long cat;            /* category */
 // long c;              /* number to encode in each category */
 // int runlength_0,runlength_1;      /* run length */
  long symbols;        /* number of symbols to encode */
  long *cache_sym;     /* temporary storage for symbols */
 // long *cache_c;       /* temporary storage for numbers */
  long pred_error;     /* prediction error for DC coefficients */

  /* allocate memory */
  alloc_long_vector(&cache_sym,(nx+2)*(ny+2));
 // alloc_long_vector(&cache_c,nx*ny);
  


  /* initialise symbol counter */
/*  symbols = 0;
  runlength_0=-1;
  runlength_1=-1;
  for(k=0;k<nx;k++)
  {
    for(l=0;l<ny;l++)
    {
      if(image[k][l]==0)
      {
        if(runlength_1>1) 
          fprintf(debug_file," 1%d", runlength_1);
        else
          for(i=runlength_1;i>-1;i--)
            fprintf(debug_file," 1");
        runlength_0++;
        runlength_1=-1;
        
      }
      else if(image[k][l]>0)
        {
          if(runlength_0>1)
            if(debug_file!=0)
          fprintf(debug_file," 0%d", runlength_0);
           else
            for(i=runlength_0;i>-1;i--)
              if(debug_file!=0)
              fprintf(debug_file," 0");
          runlength_1++;
          runlength_0=-1;

      }
  }
  }*/
  printf("\n here\n");
int flag=0,num_sym=0,run=0;
char c;
   for(k=1;k<=nx;k++)
  {
    for(l=1;l<=ny;l++)
    {
      if(image[k][l]==0)
      {
        run++;
      }
      else
      {
       cache_sym[i]=run;
       if(run>num_sym)
        num_sym=run;

     //  printf("%d %ld ", i,cache_sym[i]);
    /*   flag=0;
       for(j=0;j<i;j++)
       {
          if(cache_sym[j]==cache_sym[i])
            flag=1;

       }*/
       i++;

     /*  if(flag==0)
        num_sym++;*/
       run=0;



       
      }
  }
}
//printf(" \n %ld %d ", i, num_sym);
write_long_bitwise(nx,sizeof(nx)*8,debug_file,fp);
write_long_bitwise(ny,sizeof(ny)*8,debug_file,fp);
printf("%ld ", sizeof(long) );
encode_adaptive_wnc(cache_sym,    i,             num_sym+1,             0.3,           8192,              0,   fp);
  /* free memory */

  disalloc_long_vector(cache_sym,nx*ny);
//  disalloc_long_vector(cache_c,nx*ny);

  return;
  
}
/*--------------------------------------------------------------------------*/
void RLE_decode(long  **image,
                  FILE* debug_file,
                  BFILE* fp)
{
  long nx,ny;
  read_long_bitwise(&nx,sizeof(long)*8,debug_file,fp);
  read_long_bitwise(&ny,sizeof(long)*8,debug_file,fp);
  printf("%ld %ld \n",nx,ny );
}

/*--------------------------------------------------------------------------*/
float mse(long **u, long **f, long nx, long ny) {
  long i,j,c;
  long sum, diff;

  sum=0;

    for (i=1;i<=nx;i++)
      for (j=1;j<=ny;j++) {
        diff = u[i][j]-f[i][j];
        sum += diff*diff;
      }
  
  return (double)sum/(double)(nx*ny);
}





/*--------------------------------------------------------------------------*/
/*---------------- TREE STRUCTURE ------------------------------------------*/

/* definition of datatype reednode */
typedef struct REEDNode REEDNode;

/* definition of datatype reedtree */
typedef struct REEDTree REEDTree;

/* definition of a R-EED Node */
struct REEDNode {
    long value;                      /* 0: leaf, 1: internal,
                                      -1: uninitialised */
    float error;                    /* cached error for inpainting */
    long run[11];
    long  cx,cy;                    /* corner locations of the corresponding
                                       subimage */
    long  nx,ny;                    /* dimensions of the corresponding
                                       subimage */
    struct REEDNode *parent;        /* parent node */
    struct REEDNode *children[ 2 ]; /* child nodes */
    struct REEDTree *tree;          /* tree this node belongs to */
    unsigned long level;            /* tree level the node is located at,
                                       the root node has level 0 */
};

/* definition of a R-EED Tree */
struct REEDTree {
    unsigned long level;        /* number of level the tree has */
    unsigned long max_level;    /* maximum level of the tree */
    unsigned long min_level;    /* minimum level of the tree */
    long *open_nodes;           /* unused nodes per level */
    struct REEDNode **nodes;    /* all tree nodes, sorted by level */
};

/*--------------------------------------------------------------------------*/
int createREEDTree
( REEDTree *tree,               /* pointer to tree */
  const unsigned long level )   /* number of desired tree level */
/* create a r-eed tree return 0 on success, 1 else */
{
  unsigned long i = 0, j = 0;     /* loop variables */
  unsigned long n = 0;            /* number of nodes (level) */

  tree->max_level = 0;
  tree->min_level = 0;

  /* validate pointer */
  if( tree == NULL )
    return 1;

  /* allocate memory for all desired levels */
  tree->nodes = (REEDNode**)malloc( sizeof( REEDNode* ) * level );
  tree->open_nodes = (long*)malloc( sizeof( long ) * level );
  if( tree->nodes == NULL )
    return 1;

  /* allocate nodes on each level */
  for( i = 0; i < level; i++ ) {

    /* estimate the desired number of nodes for the current level */
    n = pow( 2.0, i );

    tree->nodes[ i ] = (REEDNode*)malloc( sizeof( REEDNode ) * n );
    if(tree->nodes[i]==NULL)
      printf("\n NULL" );

    if( tree->nodes[ i ] == NULL ) {

      for( j = 0; j < i; j++ )
        free( tree->nodes[ i ] );

      free( tree->nodes );
      tree->nodes = NULL;

      return 1;

    }

    tree->open_nodes[i]=n;
    
    /* initialize nodes on current level */
    if( i == 0 ) {
      /* root node */
      tree->nodes[ i ][ 0 ].value = 0;
      tree->nodes[ i ][ 0 ].parent = NULL;
      tree->nodes[ i ][ 0 ].children[ 0 ] = NULL;
      tree->nodes[ i ][ 0 ].children[ 1 ] = NULL;
      tree->nodes[ i ][ 0 ].tree = tree;
      tree->nodes[ i ][ 0 ].level = 0;
    }
    else {

      for( j = 0; j < n; j++ ) {

        /* current node */
        tree->nodes[ i ][ j ].value = -1;
        tree->nodes[ i ][ j ].parent = &tree->nodes[ i - 1 ][ j / 2 ];
        tree->nodes[ i ][ j ].children[ 0 ] = NULL;
        tree->nodes[ i ][ j ].children[ 1 ] = NULL;
        tree->nodes[ i ][ j ].tree = tree;
        tree->nodes[ i ][ j ].level = i;

        /* assignment (as child) to parent node */
        tree->nodes[ i - 1 ][ j / 2 ].children[ j % 2 ] =
          &tree->nodes[ i ][ j ];
      }

    }
 //   printf("\n");
 //   for( j = 0; j < n; j++ )
 //     printf("%d ", tree->nodes[i][j].value);

  }

  /* set the number of tree level */
  tree->level = level;

  return 0;
}

/*--------------------------------------------------------------------------*/
void set_node(REEDNode *node,
              long value,
              long cx, long cy,
              long nx, long ny,
              float error) {

  node->value = value;
  //
  node->cx = cx;
  node->cy = cy;

  node->nx = nx;
  node->ny = ny;
  node->error = error;

}

/*--------------------------------------------------------------------------*/
void init_tree( REEDTree *tree, long  **u, long threshold, long nx, long ny,
FILE* debug_file, BFILE* fp, long *num_symbols, long *source_length) {

  long i,n,j,k,l,c,prev,flag=0;  /* loop variables */
  long cx,cy; 
  long num_sym=0,num_sym_1=0,num_sym_2=0,num_sym_3=0,master_count=0,master_count_1=0,master_count_2=0,master_count_3=0;
  long run_temp; /* corners of subimage */
  long nx_new=0, ny_new=0; /* dimensions of subimage */
  REEDNode node,node1,node2;     /* current node */
  REEDNode parent;     /* parent node */
  float avg_thres=0;

  long *cache_sym,*cache; 
alloc_long_vector(&cache_sym,20000);
long *cache_sym_1,*cache_sym_2,*cache_sym_3;
alloc_long_vector(&cache_sym_1,200);
alloc_long_vector(&cache_sym_2,9000);
alloc_long_vector(&cache_sym_3,5000);
//printf("%ld \n", pow(2,10));

  /* initialise root node with full image */
  set_node(&(tree->nodes[0][0]),1,1,1,nx,ny,-1);

  /*for(i=1;i<=nx;i++)
    for(j=1;j<=ny;j++)
      avg_thres+=u[i][j];
  avg_thres/=(nx*ny);
  i=0;j=0;*/
/*  while(avg_thres>threshold)
  {
    i++;
    n = pow( 2.0, i );
    for( j = 0; j < n; j++ ) 
    {
        node = tree->nodes[i][j];
        parent = *(node.parent);
    }


  }*/

  for( i = 1; i < (long)(tree->level); i++ ) {
    n = pow( 2.0, i );
      for( j = 0; j < n; j++ ) {
        
        node = tree->nodes[i][j];
        // printf("\n %d \n",);
      //  fprintf(debug_file, "2 ");
        parent = (tree->nodes[i-1][(int)floor(j/2)]);
      //  printf("%d \n", parent.value);
       // fprintf(debug_file,"%d %d %d %d %d \n",parent.cx,parent.cy,parent.nx, parent.ny,parent.value );
        if(parent.value==1)
        {
        /* split subimage in its largest dimension */
      /* left child */
         
           
            


        if ((j%2)==0) { 
          cx=parent.cx; cy=parent.cy;
          if (parent.nx>=parent.ny) {
            nx_new=parent.nx/2;
            ny_new=parent.ny;
          } else {
            nx_new = parent.nx;
            ny_new = parent.ny/2;
          }
        } 
        /* right child */
        else { 
          if (parent.nx>=parent.ny) {
            cx=parent.cx+parent.nx/2;
            nx_new=parent.nx-parent.nx/2;
            cy = parent.cy; ny = parent.ny;
          } else {
            cy=parent.cy+parent.ny/2;
            ny_new = parent.ny-parent.ny/2;
            cx=parent.cx; nx=parent.nx;
          }
        }
        //fprintf(debug_file, "\n (%ld,%ld) (%ld,%ld) (%ld,%ld) (%ld,%ld)",cx,cy,cx+nx_new-1,cy,cx,cy+ny_new-1,cx+nx_new-1,cy+ny_new-1);
        int count=0;
        for(k=cx;k<cx+nx_new;k++)
          for(l=cy;l<cy+ny_new;l++)
            if(u[k][l]>0)
              {
                avg_thres+=u[k][l];
                count++;
              }
        if(count<=10)
         {  
          set_node(&(tree->nodes[i][j]),0,cx,cy,nx_new,ny_new,-1);
        //  continue;
          run_temp=0;
          c=0;
          for(k=cx;k<cx+nx_new;k++)
          for(l=cy;l<cy+ny_new;l++)
            { if(u[k][l]==0)
                run_temp++;

              if(u[k][l]>0)
              {
               // avg_thres+=u[k][l];
              // fprintf(debug_file, " \n %ld %ld ",k,l);
                node.run[c]=run_temp;
                c++;
                run_temp=0;
              }
            }
          //  fprintf(debug_file, "\n");
            node.run[c]=-1;
          //  printf("\n %d %d %d %d",cx,cy,nx_new,ny_new );
          //  master_count+=4;
          //  fprintf(debug_file,"\n %ld %ld->",i,j);
            for(k=0;k<=c;k++)
            {
              cache_sym[master_count]=node.run[k]+1;
              if(node.run[k]>num_sym)
                num_sym=node.run[k]+1;
              //flag=0;
              //for(prev=0;prev<master_count;prev++)
              //  if(cache_sym[prev]==cache_sym[master_count])
             //     {flag=1;
              //      break;}
             //     if(flag==0)
            //        num_sym++;
              master_count++;
            //  fprintf(debug_file," %ld",node.run[k] );
              if(i<=8)
              {
              	cache_sym_1[master_count_1]=node.run[k]+1;
              	if(node.run[k]>num_sym_1)
                	num_sym_1=node.run[k]+1;
                master_count_1++;

              }
              else if(i>8 )
              {
              	cache_sym_2[master_count_2]=node.run[k]+1;
              	if(node.run[k]>num_sym_2)
                	num_sym_2=node.run[k]+1;
                master_count_2++;
              }
            /*  else
              {
              	cache_sym_3[master_count_3]=node.run[k]+1;
              	if(node.run[k]>num_sym_3)
                	num_sym_3=node.run[k]+1;
                master_count_3++;
              }*/
         }

          }
         else
            set_node(&(tree->nodes[i][j]),1,cx,cy,nx_new,ny_new,-1);
    /*    avg_thres/=(nx_new*ny_new);
        printf("%d %d count %d avg_thres  %f \n", i,j, count, avg_thres);
        if(avg_thres>threshold)
          
          {
             
            node.value=1;
            set_node(&(tree->nodes[i][j]),1,cx,cy,nx_new,ny_new,-1);
            //printf("%d \n", node.value);
        //    printf("%d %d ", node1.value, node2.value);

          }
          else
             set_node(&(tree->nodes[i][j]),0,cx,cy,nx_new,ny_new,-1);*/
       
      }

        
      }
      printf("%ld \n",master_count );
      
    }
  /*  long num_sym_2=0;
    alloc_long_vector(&cache,num_sym);
    for(k=0;k<master_count;k++)
     { 
        
        flag=0;
        for(l=0;l<num_sym_2;l++)
        {
          if(cache_sym[k]==cache[l])
            flag=1;
        }
        if(flag==0)
        {
          cache[num_sym_2]=cache_sym[k];
       //   ret_cache[num_sym_2]=cache[num_sym_2];
          num_sym_2++;
        }

      }
      cache[num_sym_2]=num_sym_2;*/
      
 //   for(int i=0;i<=num_sym;i++)
  //    cache_sym[2+i]=cache[i];
 

 /*     for(k=0;k<master_count;k++)
      {
        for(l=0;l<num_sym_2;l++)
          if(cache_sym[k]==cache[l])
           { 
          //  printf("%ld->%ld ",cache_sym[k],l);
            cache_sym[k]=l;
            break;
            }
      } */

    
     
  // for(k=0;k<master_count;k++)
  //  printf("%ld ",cache_sym[k] );  


  //  printf("here \n \n \n \n");

    

    
   //printf( "\n %ld  %ld", *source_length ,*num_symbols);
    


    for(i=(long)(tree->level)-1;i>=1;i--)
    {
      n = pow( 2.0, i );
      k=0;
      for( j = 0; j < n; j++ ) {
        node=tree->nodes[i][j];
        k+=node.value;

      }
      if(k==-n)
        {
          free(tree->nodes[i]);
          tree->level--;
        }
      else
        break;
    }
    // printf("%d ",tree->level);
     tree->max_level=tree->level;

    

     for(i=(long)(tree->level)-1;i>=1;i--)
     {
        n = pow( 2.0, i );
      k=0;
      for( j = 0; j < n; j=j+2 ) {
        node1=tree->nodes[i][j];
        node2=tree->nodes[i][j+1];
        k+=(node1.value+node2.value)%2;

      }
      if(k==0)
        tree->max_level--;
      else
        break;

     }
     tree->min_level=0;
     for(i=1;i<=(long)tree->level;i++)
     {
         n = pow( 2.0, i );
      k=0;
      for( j = 0; j < n; j=j+2 ) {
        node1=tree->nodes[i][j];
        node2=tree->nodes[i][j+1];
        k+=(node1.value+node2.value)%2;

      }
      if(k==0)
        tree->min_level++;
      else
        break;
     }
     printf("here \n \n \n \n");
     //encode_adaptive_wnc(cache_sym,    master_count,             num_sym+1,             0.3,           8192,             debug_file,   fp);
     encode_adaptive_wnc(cache_sym_1,    master_count_1,             num_sym_1+1,             0.3,           8192,             debug_file,   fp);
   encode_adaptive_wnc(cache_sym_2,    master_count_2,             num_sym_2+1,             0.3,           8192,             debug_file,   fp);
   // encode_adaptive_wnc(cache_sym_3,    master_count_3,             num_sym_3+1,             0.3,           8192,             debug_file,   fp);
    *num_symbols=num_sym+1;
    *source_length=master_count;
    // printf("%d %d %d",tree->level, tree->min_level, tree->max_level);

 /*     for(i=0;i<=(long)tree->level;i++)
       {
       printf("\n");
       n = pow( 2.0, i );
    for( j = 0; j < n; j++ )
    printf("%d %d %d %d %d %d %d \n", i,j,tree->nodes[i][j].value,tree->nodes[i][j].cx,tree->nodes[i][j].cy,tree->nodes[i][j].nx,tree->nodes[i][j].ny);
}*/
     
}

/*--------------------------------------------------------------------------*/

int destroyREEDTree

    ( REEDTree *tree )              /* pointer to tree */

    /* destroy a r-eed tree
       return 0 on success, 1 else */

{
    unsigned long i = 0;            /* loop variable */

    /* validate pointer */
    if( tree == NULL )
        return 1;

    if( tree->nodes == NULL ) {
        tree->level = 0;
        return 0;
    }

    /* destroy nodes array */
    for( i = 0; i < tree->level; i++ )
        free( tree->nodes[ i ] );
    free( tree->nodes );
    free( tree->open_nodes);

    tree->nodes = NULL;
    tree->level = 0;

    return 0;
}

/*--------------------------------------------------------------------------*/

REEDNode* addREEDNodeByPos

    ( REEDTree *tree,             /* pointer to parent node */
      long level,                 /* level */
      long pos)                   /* position in level */

    /* add a child node by given parent node,
       return 0 on success, 1 else. */

{
    REEDNode *node = NULL;          /* temporary node */
    long n;

    /* validate pointer */
    if( tree == NULL )
        return 0;

    if( (level < 0) || (level > (long)(tree->level)-2))
      return 0;
   
    /* check if there are still open nodes */
    if (tree->open_nodes[level] < 1) {
      return 0;
    }

    /* check if position is valid */
    n = pow( 2.0, level );

    if ( !(pos >= 0 || pos < n))
      return 0;

    /*
    if( !( level >= 0 || level <= tree->max_level ) ||
        !(pos >= 0 || pos < n))
        return 0;
    */

    node = &(tree->nodes[level][pos]);

    /* do nothing if node is already set */
    if (node->value==1) {
      return 0;
    }

    /* set desired child node as leaf node */
    node->value = 1;
    node->children[0]->value = 0;
    node->children[1]->value = 0;
    
    tree->open_nodes[level]--;

    /* validate parent has capacities or is not initialised yet */
    if( (node->parent != NULL) && node->parent->value != 1 ) {

        node = node->parent;

        /* automatically initialize all parents as internal nodes */
        while( node != NULL ) {

            /* set parents value to 1 */
            node->value = 1;
            tree->open_nodes[node->level]--;
            if (tree->open_nodes[node->level]==0) tree->min_level=level+1;

            /* verify children have either a value of 1 or 0 */
            if( node->children[ 0 ]->value == -1  )
                node->children[ 0 ]->value = 0;
            if( node->children[ 1 ]->value == -1  )
                node->children[ 1 ]->value = 0;

            node = node->parent;

        }

    }

    /* update tree levels */
    if ((long)tree->max_level < level) {
      tree->max_level=level;
    }
    if (tree->open_nodes[level]==0) tree->min_level=level+1;

    return node;
}


typedef struct _BITBUFFER {
  unsigned char *b;
  long max_size;
  long byte_pos;
  long bit_pos;
} BITBUFFER;

/*--------------------------------------------------------------------------*/

void alloc_bitbuffer
     (BITBUFFER* buffer,   /* bitbuffer */
      long  n)              /* size */
     /* allocates memory for a bitbuffer of size n */
{
  long i;
  buffer->b = (unsigned char *) malloc (n * sizeof(unsigned char));
  if (buffer->b == NULL)
    {
      printf("alloc_bitbuffer: not enough memory available\n");
      exit(1);
    }
  buffer->max_size=n;
  for (i=0;i<buffer->max_size;i++) {
    buffer->b[i]=0;
  }
  buffer->byte_pos=0;
  buffer->bit_pos=0;
return;
}

/*--------------------------------------------------------------------------*/

long bitbuffer_addbit(BITBUFFER* buffer, unsigned char bit) {
  if ((buffer->byte_pos < buffer->max_size) && (buffer->bit_pos < 8)) {
    if (bit) buffer->b[buffer->byte_pos] |= (1 << buffer->bit_pos);
    else buffer->b[buffer->byte_pos] &= ~(1 << buffer->bit_pos);
    (buffer->bit_pos)++;
    if (buffer->bit_pos > 7) {
      buffer->bit_pos = 0;
      buffer->byte_pos++;
    }
    return 1;
  }
  return 0; /* adding not successful, memory full */  
}

/*--------------------------------------------------------------------------*/

void bitbuffer_writefile(BITBUFFER* buffer, FILE* bitfile) {
  fwrite (buffer->b, sizeof(unsigned char), buffer->byte_pos+1, bitfile);
}

/*--------------------------------------------------------------------------*/

void bitbuffer_printbits(BITBUFFER* buffer) {
  long i,j;
  printf("Bitbuffer stores %ld bits:\n",
         (buffer->byte_pos)*8+buffer->bit_pos);
  for (i=0;i<buffer->byte_pos;i++) {
    for (j=0;j<8;j++) {
      printf("%u",((buffer->b[i] & (1u << j)) > 0));
    }
    printf(" ");
  }
  i=buffer->byte_pos;
  for (j=0;j<buffer->bit_pos;j++) {
    printf("%u",((buffer->b[i] & (1u << j)) > 0));
  }
  printf("\n");
}




/*--------------------------------------------------------------------------*/

int store_tree
    ( const REEDTree *tree,        /* pointer to tree */
      FILE *fp, long nx, long ny, long num_sym, long source_length, long *cache)                  /* bitfile  to write to */
    /* write reed tree to bitfile,
       return 0 on success, 1 else. */
{   
    unsigned long i = 0, j = 0;     /* loop variables */
    unsigned long n = 0;            /* number of nodes on current level */
    BITBUFFER buffer;


    /* validate pointer */
    if( tree == NULL || fp == NULL )
       return 1;

    if( tree->nodes == NULL && tree->level != 0 )
        return 1;

    /* 1. header */
    
      fwrite(&nx,sizeof(long),1,fp);
    fwrite(&ny,sizeof(long),1,fp);

    fwrite(&num_sym,sizeof(long),1,fp);
    fwrite(&source_length,sizeof(long),1,fp);

  //  for(int i=0;i<=num_sym;i++)
  //    fwrite(&(cache[i]),sizeof(long),1,fp);

    fwrite((unsigned char*)&(tree->min_level),sizeof(char),1,fp);
    fwrite((unsigned char*)&(tree->max_level),sizeof(char),1,fp);


    /* printf("Writing tree (%ld %ld)\n",tree->min_level,tree->max_level); */
    
    /* 2. data */
    /* write information of valid nodes to buffer array */
    alloc_bitbuffer(&buffer,(long)(pow(2.0,tree->max_level)/8)+1);
    for( i = tree->min_level; i <= tree->max_level; i++ ) {
        /* estimate number of nodes on current level */
        n = pow( 2.0, i );
        for( j = 0; j < n; j++ ) {
            if( (tree->nodes[i][j].value) != -1 ) {
              
              /*printf("wrote bit %ld for node %ld %ld\n",
                tree->nodes[ i ][ j ].value,i,j);*/
              
              /*bfputb( tree->nodes[ i ][ j ].value, fp );*/
              bitbuffer_addbit(&buffer,
                               (unsigned char)(tree->nodes[i][j].value));
            }
        }
    }



    /*    bitbuffer_printbits(&buffer); */
   bitbuffer_writefile(&buffer,fp);



    return 0;
}

/*--------------------------------------------------------------------------*/

int load_tree
    ( REEDTree *tree,        /* pointer to tree */
      FILE *fp, BFILE *fp2, FILE *debug_file)              /* bitfile  to write to */
    /* write reed tree to bitfile,
       return 0 on success, 1 else. */
{
    unsigned long l,k,i = 0, j = 0,master_count=0;     /* loop variables */
    unsigned long n = 0;            /* number of nodes on current level */
    long bit = 0;  
    long num_sym,source_length;                  /* single bit, being loaded from bitfile */
    unsigned char byte;
    long bit_pos;
    REEDNode node,parent;
    long nx,ny,cx,cy,nx_new,ny_new;
    long *cache_sym,*cache;


    /* validate pointer */
    if( tree == NULL || fp == NULL )
       return 1;

    /* 1. header */
     fread(&nx,sizeof(long),1,fp);
    fread(&ny,sizeof(long),1,fp);
    fread(&num_sym,sizeof(long),1,fp);
    fread(&source_length,sizeof(long),1,fp);
    alloc_long_vector(&cache,num_sym+1);
  //  for(int i=0;i<=num_sym;i++)
  //    fread(&(cache[i]),sizeof(long),1,fp);
    //printf("\n %ld %ld %ld %ld", nx,ny, num_sym, source_length);

    fread (&byte, sizeof(char), 1, fp);
    tree->min_level = (long)byte;
    fread (&byte, sizeof(char), 1, fp);
    tree->max_level = (long)byte;
    alloc_long_vector(&cache_sym,source_length);
    decode_adaptive_wnc(fp2, source_length, num_sym,0.3,8192,0, cache_sym);
    cache_sym[source_length-1]=0;


    for(k=0;k<source_length;k++)
    {
      cache_sym[k]=cache_sym[k]-1;    
     
    }
   // for(i=0;i<source_length;i++)
   //   printf("%ld ", cache_sym[i]);

    

    /*printf("loading tree (min %ld, max %ld)\n",tree->min_level,
      tree->max_level);*/

    /* 2. data */

    /* split all nodes below minimum level */
    for( i = 0; i < tree->min_level; i++ ) {
        /* estimate number of nodes on current level */
        n = pow( 2.0, i );

        tree->open_nodes[i]=0;
        for( j = 0; j < n; j++ ) {
            tree->nodes[i][j].value=1;
        }
    }

    bit_pos = 8;
    /* load all nodes up to maximum level */
    for( i = tree->min_level; i <= tree->max_level; i++ ) {
            /* estimate number of nodes on current level */
            n = pow( 2.0, i );

            for( j = 0; j < n; j++ ) {
              /* check for each node if it can exist (i.e. parent is split) */
              if (tree->nodes[i][j].parent->value == 1) {
                if (bit_pos > 7) {
                  fread (&byte, sizeof(char), 1, fp);
                  /*printf("read byte %u\n",byte);*/
                  bit_pos = 0;
                }
                bit = (long)((byte & (1u << bit_pos)) > 0);
                /* printf("read bit %ld for node %ld %ld (bitpos %ld)\n",bit,i,j,
                   bit_pos); */
                bit_pos++;
                tree->nodes[i][j].value=bit;
                if (bit == 1) {
                  tree->open_nodes[i]--;
             }
           }
        }
     }

    /* set children on maximum level */
   i = tree->max_level + 1;
   /* estimate number of nodes on current level */
   n = pow( 2.0, i );

   for( j = 0; j < n; j++ ) {
   /* check for each node if parent is split */
     if (tree->nodes[i][j].parent->value == 1) {
       tree->nodes[i][j].value=0;
     }
   }
   // shaving off unwanted levels
   for(i=(long)(tree->level)-1;i>=1;i--)
    {
      n = pow( 2.0, i );
      k=0;
      for( j = 0; j < n; j++ ) {
        node=tree->nodes[i][j];
        k+=node.value;

      }
      if(k==-n)
        {
          free(tree->nodes[i]);
          tree->level--;
        }
      else
        break;
    }

    set_node(&(tree->nodes[0][0]),tree->nodes[0][0].value,1,1,nx-2,ny-2,-1);
for( i = 1; i < (long)(tree->level); i++ ) {
    n = pow( 2.0, i );
      for( j = 0; j < n; j++ ) {
        
        node = tree->nodes[i][j];
        // printf("\n %d \n",);
      //  fprintf(debug_file, "2 ");
        parent = (tree->nodes[i-1][(int)floor(j/2)]);
      //  printf("%d \n", parent.value);
       // fprintf(debug_file,"%d %d %d %d %d \n",parent.cx,parent.cy,parent.nx, parent.ny,parent.value );
        if(parent.value==1)
        {
        /* split subimage in its largest dimension */
      /* left child */
         
           
            


        if ((j%2)==0) { 
          cx=parent.cx; cy=parent.cy;
          if (parent.nx>=parent.ny) {
            nx_new=parent.nx/2;
            ny_new=parent.ny;
          } else {
            nx_new = parent.nx;
            ny_new = parent.ny/2;
          }
        } 
        /* right child */
        else { 
          if (parent.nx>=parent.ny) {
            cx=parent.cx+parent.nx/2;
            nx_new=parent.nx-parent.nx/2;
            cy = parent.cy; ny = parent.ny;
          } else {
            cy=parent.cy+parent.ny/2;
            ny_new = parent.ny-parent.ny/2;
            cx=parent.cx; nx=parent.nx;
          }
        }
      }
      //if(node.value==0)
      //printf("\n %ld %ld %ld %ld value= %ld->", node.cx,node.cy,node.nx,node.ny, node.value);
      set_node(&(tree->nodes[i][j]),node.value,cx,cy,nx_new,ny_new,-1);
      node = tree->nodes[i][j];
      //if(node.value==0)
      //if(i<5)
     // printf("\n %ld %ld %ld %ld %ld %ld %ld %ld", node.cx,node.cy,node.nx,node.ny, node.value,parent.cx,parent.cy,parent.nx,parent.ny);
      //printf(" %ld", node.value);
      if(node.value==0)
      {
        
        k=0;
        while(cache_sym[master_count]!=-1 && master_count<source_length)
        {
          (tree->nodes[i][j]).run[k]=cache_sym[master_count];
         // printf(" %ld", node.run[k]);
          master_count++;
          k++;
        }
        if(cache_sym[master_count]==-1)
          {
            (tree->nodes[i][j]).run[k]=-1;
            //
            master_count++;
            k++;
          }
        // printf("\n %ld %ld->", i,j);
         // for(k=0;k<11;k++)
         //   printf(" %ld",(tree->nodes[i][j]).run[k] );
      }
      
      //if(node.value==0)
     //   printf("\n %d %d %d %d",cx,cy,nx_new,ny_new );

    }
  }

 /* for( i = 0; i < (long)(tree->level); i++ ) {
    n = pow( 2.0, i );
      for( j = 0; j < n; j++ ) {

        node = tree->nodes[i][j];
        if(node.value==0)
        {
          //k=0;
         //  printf("\n %ld %ld %ld %ld", node.cx,node.cy,node.nx,node.ny);
           printf("\n %ld %ld->", i,j);
          for(k=0;k<11;k++)
            printf(" %ld",node.run[k] );
        //  while(node.run[k]!=-1)
        //    {printf(" %ld",node.run[k] );k++;}
        //  printf(" %ld", node.run[k]);
        }
      }
    } */


    return 0;
}






void image_reconstruct(long **u, REEDTree *tree)
{
  long i,j,k,n,x,y,csum=0;
  REEDNode node;
  
  for(i=0;i<tree->level;i++)
  {
  //i=8;
    n = pow( 2.0, i );
    //printf("level %ld -> %ld elements \n", i,(sizeof(tree->nodes[i])));;
      for( j = 0; j < n; j++ ) 
      {
        
      //  printf("\n %ld %ld ",i,j);
      //  printf("here 1");
        node=tree->nodes[i][j];

        if(node.value==-1 || node.value==1)
          continue;

         if(node.value==0)
        {
          k=0;
          csum=0;
          
          while(node.run[k]!=-1)
          {
           
            csum+=node.run[k];
            x=node.cx+(long)(csum/node.ny);
            y=node.cy+(long)(csum%node.ny);
           // printf("(%ld %ld (%ld, %ld), ",csum,node.run[k],x,y);
            u[x][y]=255;
            k++;
            csum++;
          }
        }

       // printf("here 2");
      }
  }
}

/*--------------------------------------------------------------------------*/
/*-------- TREE SPECIFIC FUNCTIONS -----------------------------------------*/

void alloc_vector_nodes
     (REEDNode ***vector,   /* vector */
      long  n1)         /* size */
     /* allocates memory for a vector of size n1 */
{
  *vector = (REEDNode**) malloc (n1 * sizeof(REEDNode*));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough memory available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/
void add_points(REEDNode *node, float **mask) {
  mask[node->cx][node->cy] = 1.0;
  mask[node->cx][node->cy+node->ny-1]=1.0;
  mask[node->cx+node->nx-1][node->cy]=1.0;
  mask[node->cx+node->nx-1][node->cy+node->ny-1]=1.0;
  mask[node->cx+node->nx/2][node->cy+node->ny/2]=1.0;
}


/*--------------------------------------------------------------------------*/

int main (int argc, char **args)

{
  /* user interaction */
  char   ch;                   /* for reading single chars */
  char   used[256];            /* flag for all possible chars, used to check for
                                  multiple occurences of input parameters */
  long len;
  
  /* parameters */

  /* file IO */
  char *output_file = 0;       /* file name of output image */
  char *input_file = 0;        /* file name of uncompressed image */
  char tmp_file[1000];
  char tmp_file2[1000];
  char tmp_file_run[1000];         /* string for intermediate filenames */
  char total_file[1000];       /* file name of total compressed output */
  char comments[10000];        /* string for comments */
  char *program_call;          /* call of the compression program */
  char extension[5];           /* extension of input file */
  long format;                 /* format of input file */

  /* image information */
  long nx[3], ny[3], nc;                /* image dimensions */

  /* loop variables */
  long i,j;

  /* image struct that contains all information
   * about original image, compressed image, and parameters */
  ImageData image;

  /* compression/decompression */
  long flag_compress=0;
  FILE *binary_file=0;
  BFILE *binary_file_run=0;
  BFILE *binary_file2=0;
  FILE *binary_file_in=0;
  BFILE *binary_file_in2=0;       /* binary file for compressed bitstring */
  char*  debug_file=0;        /* filename for writing debug information */
  FILE*  dfile=0;             /* file for writing debug information */
  long   q=0;                 /* quantisation parameter */
  long   s=0;                 /* chroma subsampling factor */
  long **tmp_img,**image_out,**image_out_2; 
  long *cache;            /* temporary image */
  
  printf ("\n");
  printf ("PROGRAMMING EXERCISE FOR IMAGE COMPRESSION\n\n");
  printf ("**************************************************\n\n");
  printf ("    Copyright 2017 by Pascal Peter                \n");
  printf ("    Dept. of Mathematics and Computer Science     \n");
  printf ("    Saarland University, Saarbruecken, Germany    \n\n");
  printf ("    All rights reserved. Unauthorised usage,      \n");
  printf ("    copying, hiring, and selling prohibited.      \n\n");
  printf ("    Send bug reports to                           \n");
  printf ("    peter@mia.uni-saarland.de                     \n\n");
  printf ("**************************************************\n\n");

  /* ARGUMENT PROCESSING ******************************************************/
  init_image(&image); /* initialise image with standard parameters */

  for (i = 0; i <= 255; i++) {
    used[i] = 0;
  }

  while ((ch = getopt(argc,args,"i:q:o:D:s:")) != -1) {
    used[(long)ch]++;
    if (used[(long)ch] > 1) {
      printf("Duplicate parameter: %c\n",ch);
      printf("Please check your input again.\n");
    }
    
    switch(ch) {
    case 's': s=atoi(optarg); image.s=s; break;
    case 'q': q=atoi(optarg);break;
    case 'i': input_file = optarg;break;
    case 'o': output_file = optarg;break;
    case 'D': debug_file = optarg;break;
    default:
      printf("Unknown argument.\n");
      print_usage_message();
      return 0;
    }
  }

  /* reset quantisation matrix */
  if (q>0) {
    for (i=0;i<8;i++)
      for (j=0;j<8;j++)
        w[i][j]=q;
  }

  if (s==0) s = image.s;
  
  if (output_file == 0 || input_file == 0) {
    printf("ERROR: Missing mandatory parameter, aborting.\n");
    print_usage_message();
    return 0;
  }

  /* prepare file names */
  sprintf(total_file,"%s.coded",output_file);

  /* create reboot string */
  len = 0;
  for (i=0;i<argc;i++) {
    len += strlen(args[i]) + 1;
  }

  program_call = (char*)malloc( sizeof(char) * ( len + 3 ) );
  sprintf(program_call," ");
  for (i=0;i<argc;i++) {
    sprintf(program_call, "%s%s ",program_call,args[i]);
  }

  /* DETERMINE COMPRESSION/DECOMPRESSION MODE**********************************/

  /* try to identify the file format via the extension */
  strcpy(extension, input_file+(strlen(input_file)-4));
  if      (!strcmp(extension, ".pgm")) format = FORMAT_PGM;
  else if (!strcmp(extension, ".ppm")) format = FORMAT_PPM;
  else {
    printf("ERROR: Extension %s not supported for input file, aborting.\n",
           extension);
    print_usage_message();
    return 0;
  }

  /* determine if we are in compression or decompression mode */
  if (format==FORMAT_PPM || format==FORMAT_PGM) {
    flag_compress = 1;
  } else {
    flag_compress = 0;
  }
  
  if (flag_compress == 1) {
    /* COMPRESS ***************************************************************/

    /* read input image */
    if (format==FORMAT_PPM) {
      read_ppm_and_allocate_memory(input_file,&image.nx,&image.ny,
                                   &image.orig_rgb);
      image.nc = 3;
      printf("Image %s loaded (PPM).\n\n",input_file);
    } else if (format==FORMAT_PGM) {
      read_pgm_header(input_file,&image.nx,&image.ny);
      alloc_long_cubix(&image.orig_rgb,MAXCHANNELS,image.nx+2,image.ny+2);
      read_pgm_and_allocate_memory(input_file,&image.nx,&image.ny,
                                   &image.orig_rgb[0]);
      image.nc = 1;
      printf("Image %s loaded (PGM).\n\n",input_file);
    }
    nx[0] = image.nx; ny[0] = image.ny; nc = image.nc;
    image.size_orig=get_size_of_file(input_file);

    printf("Image dimensions: %ld x %ld x %ld\n",nx[0],ny[0],nc);


    /* allocate memory */
    alloc_image(&image,nx[0],ny[0]);
    //alloc_long_matrix(&tmp_img,nx[0],ny[0]);
    
   
   // printf("\n %ld %ld",sizeof(image.orig_rgb[0][]),image.ny_ext[0]+2);
    alloc_long_matrix(&image_out,nx[0]+2,ny[0]+2);
    alloc_long_matrix(&image_out_2,nx[0]+2,ny[0]+2);

    //alloc_long_vector(&cache,nx[0]*ny[0]);



    /* convert to YCbCr space or copy over grey value image */
   /* if (nc > 1) {
      RGB_to_YCbCr(image.orig_rgb,image.orig_ycbcr,nx[0],ny[0]);
    } else {
      copy_matrix_long(image.orig_rgb[0],image.orig_ycbcr[0],nx[0],ny[0]);
    }*/
    
    /* perform chroma subsampling */
    /*
    if (nc > 1) {
      nx[1]=nx[2]=nx[0]/s;
      ny[1]=ny[2]=ny[0]/s;    
      if ((nx[0] % s) > 0) {nx[1]++;nx[2]++;}
      if ((ny[0] % s) > 0) {ny[1]++;ny[2]++;}
      if (s > 1) {
       for(i=0;i<nx[1];i++)
        {
          printf("\n");
          for(j=0;j<ny[1];j++)
          printf(" %ld",image.orig_ycbcr[1][i][j]);
        }
        printf("\n");
        for(i=0;i<nx[2];i++)
        {
          printf("\n");
          for(j=0;j<ny[2];j++)
          printf(" %ld",image.orig_ycbcr[2][i][j]);
        }*/
   /*     subsample(image.orig_ycbcr[1],image.rec_quant[0],nx[0],ny[0],s);
        copy_matrix_long(image.rec_quant[0],image.orig_ycbcr[1],image.nx_ext[1],
                         image.ny_ext[1]);
        subsample(image.orig_ycbcr[2],image.rec_quant[0],nx[0],ny[0],s);
        copy_matrix_long(image.rec_quant[0],image.orig_ycbcr[2],image.nx_ext[1],
                         image.ny_ext[1]);
      }
      printf("Chroma subsampling by factor %ld (%ld x %ld -> %ld x %ld)\n",
             s,nx[0],ny[0],nx[1],ny[1]); */
    }
    
    /* extend image dimensions to multiples of block_size */
   /* for (i=0; i<nc; i++) {
      extend_image(image.orig_ycbcr[i],nx[i],ny[i],image.block_size,
                   image.orig_ycbcr[i]);
    }*/

    /* open debug file if debug mode is active */
    if (debug_file != 0) {
      dfile = fopen(debug_file,"w");
    }

    /* prepare files for encoding */
   /* printf("Computing DCT and quantising coefficients\n");
    sprintf(tmp_file,"%s.wnc",output_file);
    binary_file = bfopen(tmp_file,"w"); */

    /* apply block DCT and encode */
  /*  for (i=0; i<nc; i++) {
      block_DCT(image.orig_ycbcr[i],image.nx_ext[i],image.ny_ext[i],
                image.block_size,image.dct[i]);
      block_quantise(image.dct[i],image.nx_ext[i],image.ny_ext[i],0,
                     image.dct_quant[i]);

      block_encode(image.dct_quant[i],image.nx_ext[i],image.ny_ext[i],
                   dfile,binary_file);          
    } */
    //tmp_file2=tmp_file; 


  /*  sprintf(tmp_file_run,"%s.run",output_file);
   binary_file = bfopen(tmp_file_run,"w");



    RLE_encode(image.orig_rgb[0],nx[0],ny[0],dfile,binary_file);
    //RLE_encode(tmp_img,nx[0],ny[0],dfile,binary_file);




    //fprintf(dfile,"\n here");
    bfclose(binary_file);



    binary_file = bfopen(tmp_file_run,"r");
    RLE_decode(image_out_2,dfile, binary_file);
    bfclose(binary_file);   */


    REEDTree tree,tree2;          //  <-------------------------------------uncomment here
    
    createREEDTree(&tree,20);
    createREEDTree(&tree2,20);

   // REEDTree *tree, long  **u, long threshold, long nx, long ny,
//FILE* debug_file
       sprintf(tmp_file,"%s.tree",output_file);
    sprintf(tmp_file2,"%s.treerun",output_file);

    binary_file = fopen(tmp_file,"w");
    binary_file2=bfopen(tmp_file2,"w");
   
   long num_sym,n;
    init_tree(&tree, image.orig_rgb[0],2,nx[0],ny[0],dfile,binary_file2, &num_sym, &n);
    store_tree(&tree, binary_file,nx[0]+2,ny[0]+2, num_sym, n, cache);

    
  /*  long a[5]={0,1,2,3,4};
    encode_adaptive_wnc(a,5,6,0.3,128,dfile,binary_file2);*/
    
    fclose(binary_file);
    bfclose(binary_file2);
   

   /*binary_file_in=fopen(tmp_file,"r");
    binary_file_in2=bfopen(tmp_file2,"r");
    load_tree(&tree2,binary_file_in,binary_file_in2, dfile);
    image_reconstruct(image_out,&tree2);


    sprintf(tmp_file,"%s.pgm",output_file);
    normalise_to_8bit(image_out,nx[0],ny[0],image_out);
    write_pgm(image_out, nx[0], ny[0],tmp_file, comments);
    printf("Resulting MSE: %f\n",mse(image.orig_rgb[0],image_out,nx[0],ny[0] ));
      

   // decode_adaptive_wnc(binary_file_in2,5,6,0.3,1024,dfile,a);
      fclose(binary_file_in);
    bfclose(binary_file_in2); */
        
 /*   long max=0;
    for(i=0;i<nx[0];i++)
      { 
        for(j=0;j<ny[0];j++)
        {
          if(image_out[i][j]>max)
            max=image_out[i][j];
        }
        
    }
    printf(" max=%ld",max );*/
    
    




    //decode_tree(&tree2, binary_file_in2 )
   // printf("%d %d %d",tree2.level, tree2.min_level, tree2.max_level);
  /*  int diff=0;
    long n;
    for( i = 0; i < tree.max_level; i++ ) {
        
        n = pow( 2.0, i );

      
        for( j = 0; j < n; j++ ) {
           diff+=( tree.nodes[i][j].value-tree2.nodes[i][j].value);
        }
    }

    printf("\n %d", diff);*/

    //printf("\n here \n");

    /* close binary file */
   // 

    /* output image information */
   /* printf("Resulting compression ratio: %f:1\n\n", 
           get_compression_ratio(input_file,tmp_file));*/
    
    /* write image data */
  /*  write_comment_string(&image,0,comments);
    for (i=0; i<nc; i++) {
      sprintf(tmp_file,"%s_dct_channel%ld.pgm",output_file,i);
      abs_img(image.dct_quant[i],nx[i],ny[i],tmp_img);
      normalise_to_8bit(tmp_img,nx[i],ny[i],tmp_img);
      write_pgm(tmp_img, nx[i], ny[i],tmp_file, comments);
    } */

    /* reconstruct */
 /*   printf("Requantising and compute inverse DCT\n");
    for (i=0; i<nc; i++) {
      block_requantise(image.dct_quant[i],image.nx_ext[i],image.ny_ext[i],0,
                     image.dct_quant[i]);
      block_IDCT(image.dct_quant[i],image.nx_ext[i],image.ny_ext[i],
                 image.block_size,image.rec[i]);
      convert_matrix_int(image.rec[i],image.rec_quant[i],
                        image.nx_ext[i],image.ny_ext[i]);
    }*/

    /* perform upsampling if downsampling was applied before */
  /*  if (s>1 && nc > 1) {
      upsample(image.rec_quant[1],tmp_img,nx[0],ny[0],s);

      copy_matrix_long(tmp_img,image.rec_quant[1],nx[0],ny[0]);
      upsample(image.rec_quant[2],tmp_img,nx[0],ny[0],s);
      copy_matrix_long(tmp_img,image.rec_quant[2],nx[0],ny[0]);

      printf("Chroma upsampling by factor %ld (%ld x %ld -> %ld x %ld)\n",
             s,nx[1],ny[1],nx[0],ny[0]);
    } */

    /* convert back from YCbCr to RGB */
  /*  if (nc>1) {
      YCbCr_to_RGB(image.rec_quant,image.rec_quant,image.nx_ext[0],
                   image.ny_ext[0]);
    }

    printf("Resulting MSE: %f\n",mse(image.orig_rgb,image.rec_quant,nx[0],ny[0],
                                     nc)); */

    /* write reconstruction */
  /*  if (format==FORMAT_PPM) {
      sprintf(tmp_file,"%s_rec.ppm",output_file);
      write_ppm(image.orig_rgb, nx[0], ny[0], tmp_file, comments);
    } else {
      sprintf(tmp_file,"%s_rec.pgm",output_file);
      write_pgm(image.orig_rgb[0], nx[0], ny[0], tmp_file, comments);
    } */
    
    if (debug_file !=0) {
      fclose(dfile);
     }
     
   //else {
    /* DECOMPRESS *************************************************************/

 // }
  
  /* ---- free memory  ---- */
  //disalloc_long_matrix(tmp_img,image.nx_ext[0]+2,image.ny_ext[0]+2);
  disalloc_long_matrix(image_out,nx[0],ny[0]);
  destroy_image(&image);
  free(program_call);

  return(0);
}


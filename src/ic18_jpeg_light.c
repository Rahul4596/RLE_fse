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
  double r,           /* rescaling parameter */
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

  /* last step */
  if (debug_file != 0) {
    fprintf(debug_file,"last interval - written bits:");
  }
  if (L<M14) {
      set_bit(compressed,0);
      if (debug_file != 0) {
        fprintf(debug_file,"0");
      }
    for (j=0;j<k+1;j++) {
      set_bit(compressed,1);
      if (debug_file != 0) {
        fprintf(debug_file,"1");
      }
    }
  } else {
    set_bit(compressed,1);
    if (debug_file != 0) {
      fprintf(debug_file,"1");
    }
    for (j=0;j<k+1;j++) {
      set_bit(compressed,0);
      if (debug_file != 0) {
        fprintf(debug_file,"0");
      }
    }
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
              "M=%ld is too small (C=%ld), setting M to %ld\n",C,M,8*C);
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
/* calulates block DCT of input image/channel */
void block_DCT(long  **f,            /* input image */
               long nx, long ny,     /* image dimensions */
               long N,               /* block size */
               double  **dct) {      /* output DCT coefficients */        
  long   x,y,u,v,k,l;   /* loop variables */
  double  pi;            /* variable pi */
  /* double  **f;           /\* extended image *\/ */
  /* double  **c;           /\* extended coefficient matrix *\/ */
  double  **basis;       /* time saver */
  double  *alpha;        /* array for coefficients */
  /* long   N=8;            /\* block size *\/ */
  long   blocks_x;       /* number of blocks in each direction */
  long   blocks_y;
  long   ox,oy;          /* block offsets */
  
  /* ---- set time savers ---- */
  pi = 2.0 * asinf (1.0);
  
  /* determine number of blocks */
  blocks_x = nx/N;
  if ((nx % N) > 0) blocks_x++;
  blocks_y = ny/N;
  if ((ny % N) > 0) blocks_y++;

  /* ---- allocate memory ---- */
  alloc_vector (&alpha, N);
  alloc_matrix (&basis, N,N);
 
  /* precompute scaling functions */
  alpha[0]=sqrt(1.0/(double)N);
  for (x=1;x<N;x++) {
    alpha[x]=sqrt(2.0/(double)N); 
  }

  /* precompute cosine basis functions */
  for (u=0; u<N; u++)
    for (x=0; x<N; x++) {
      basis[u][x] = cosf(pi/(double)N*((double)x+0.5)*(double)u);
    }

  /* Iterate over all blocks */
  for (k=0;k<blocks_x;k++)
    for (l=0;l<blocks_y;l++) {
      ox = k*N+1; oy=l*N+1; /* define block offsets */
 
      /* 2-D block DCT */
      /* supplement code here */
      /* HINT: For the block k/l, the offsets ox and oy help you to find
         the pixel coordinates in the full image. You can assume that
         each block has full size NxN, the image is extended appropriately
         elsewhere. Make use of the predefined cosine-basis and scaling
         functions. */
    	for (u=0; u<N; u++)
            for (v=0; v<N; v++)
            {
            	dct[ox+u][oy+v]=0;
            	for (x=0; x<N; x++)
        			for (y=0; y<N; y++)
        			{
        				dct[ox+u][oy+v] += alpha[u]*alpha[v]*(float)f[ox+x][oy+y]*
                               basis[u][x]*basis[v][y];
        				

        			}
            }
    	

    	

    }
 
  /* ---- free memory ---- */
  disalloc_vector (alpha, N);
  disalloc_matrix (basis, N, N);
  
  return;
}

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
/* inverts block DCT of quantised input coefficients */
void block_IDCT(long   **dct,       /* quantised input coefficients */
                long   nx, long ny,  /* image dimensions */
                long   N,            /* block_size */
                double **f) {       /* output image */        
  long   x,y,u,v,k,l;   /* loop variables */
  double  pi;            /* variable pi */
  double  **basis;       /* time saver */
  double  *alpha;        /* array for coefficients */
  long   blocks_x;       /* number of blocks in each direction */
  long   blocks_y;
  long   ox,oy;          /* block offsets */

  /* ---- set time savers ---- */
  pi = 2.0 * asinf (1.0);
  
  /* determine number of blocks */
  blocks_x = nx/N;
  if ((nx % N) > 0) blocks_x++;
  blocks_y = ny/N;
  if ((ny % N) > 0) blocks_y++;

  /* ---- allocate memory ---- */
  alloc_vector (&alpha, N);
  alloc_matrix (&basis, N,N);

  /* precompute scaling functions */
  alpha[0]=sqrt(1.0/(double)N);
  for (x=1;x<N;x++) {
    alpha[x]=sqrt(2.0/(double)N); 
  }

  /* precompute cosine basis functions */
  for (u=0; u<N; u++)
    for (x=0; x<N; x++) {
      basis[u][x] = cosf(pi/(double)N*((double)x+0.5)*(double)u);
    }

  /* Iterate over all blocks */
  for (k=0;k<blocks_x;k++)
    for (l=0;l<blocks_y;l++) {
      ox = k*N+1; oy=l*N+1; /* define block offsets */
 
      /* 2-D block DCT */
      for (x=0; x<N; x++)
        for (y=0; y<N; y++) {
          f[ox+x][oy+y]=0;
          for (u=0; u<N; u++)
            for (v=0; v<N; v++) {
              f[ox+x][oy+y] += alpha[u]*alpha[v]*dct[ox+u][oy+v]*
                               basis[u][x]*basis[v][y];
            }
        }
    }
  
  /* ---- free memory ---- */
  disalloc_vector (alpha, N);
  disalloc_matrix (basis, N, N);
  
  return;
}


/*--------------------------------------------------------------------------*/
/* quantises DCT coefficients in blocks */
void block_quantise(double  **dct,      /* input DCT coefficients */
                    long nx, long ny,   /* image dimensions */
                    FILE* debug_file,   /* 0 - no output, 
                                           otherwise debug output to file */
                    long  **quant) {    /* output quantised DCT coefficients */

  long   u,v,k,l;   /* loop variables */
  long   N=8;            /* block size */
  long   blocks_x;       /* number of blocks in each direction */
  long   blocks_y;
  long   ox,oy;          /* block offsets */
  long max, min;
  max = 0;
  min = 10000000000;
  
  /* determine number of blocks */
  blocks_x = nx/N;
  if ((nx % N) > 0) blocks_x++;
  blocks_y = ny/N;
  if ((ny % N) > 0) blocks_y++;

  if (debug_file != 0) {
    fprintf(debug_file,"QUANTISATION\n");
    fprintf(debug_file,"============\n");
  }

  /* Iterate over all blocks and apply quantisation matrix */
  for (k=0;k<blocks_x;k++)
    for (l=0;l<blocks_y;l++) {
      ox = k*N+1; oy=l*N+1; /* define block offsets */

      if (debug_file != 0) {
        fprintf(debug_file,"Block %ld %ld\n",k,l);
        fprintf(debug_file,"-------------\n");
      }
      
      /* apply quantisation matrix (handle incomplete blocks correctly!) */
      for (u=0; u<N; u++)
        for (v=0; v<N; v++) {
          if ((ox+u<=nx) && (oy+v<=ny)) {
            quant[ox+u][oy+v] = (long)(dct[ox+u][oy+v]/(double)w[u][v]);
          }
          if (debug_file != 0) {
            fprintf(debug_file,"%f -> %ld (w %ld)\n",dct[ox+u][oy+v],
                    quant[ox+u][oy+v],w[u][v]);
          }
          if (quant[ox+u][oy+v]>max) max = quant[ox+u][oy+v];
          if (quant[ox+u][oy+v]<min) min = quant[ox+u][oy+v];
        }
    }

  if (debug_file != 0) {
    fprintf(debug_file,"quantised coefficients are in range [%ld,%ld]\n",
           min,max);
  }
  
  return;
  
}

/*--------------------------------------------------------------------------*/
/* quantises DCT coefficients in blocks with INVERSE of quantisation weights */
void block_requantise(long  **quant,      /* input DCT coefficients */
                      long nx, long ny,   /* image dimensions */
                      FILE* debug_file,   /* 0 - no output, 
                                             otherwise debug output to file */
                      long  **dct) {    /* output quantised DCT coefficients */
  
  long   u,v,k,l;   /* loop variables */
  long   N=8;            /* block size */
  long   blocks_x;       /* number of blocks in each direction */
  long   blocks_y;
  long   ox,oy;          /* block offsets */
  long max, min;
  max = 0;
  min = 10000000000;
  
  /* determine number of blocks */
  blocks_x = nx/N;
  if ((nx % N) > 0) blocks_x++;
  blocks_y = ny/N;
  if ((ny % N) > 0) blocks_y++;

  if (debug_file != 0) {
    fprintf(debug_file,"QUANTISATION\n");
    fprintf(debug_file,"============\n");
  }

  /* Iterate over all blocks and apply quantisation matrix */
  for (k=0;k<blocks_x;k++)
    for (l=0;l<blocks_y;l++) {
      ox = k*N+1; oy=l*N+1; /* define block offsets */

      if (debug_file != 0) {
        fprintf(debug_file,"Block %ld %ld\n",k,l);
        fprintf(debug_file,"-------------\n");
      }
      
      /* apply quantisation matrix (handle incomplete blocks correctly!) */
      for (u=0; u<N; u++)
        for (v=0; v<N; v++) {
          if ((ox+u<=nx) && (oy+v<=ny)) {
            dct[ox+u][oy+v] = (long)(quant[ox+u][oy+v]*(double)w[u][v]);
          }
          if (debug_file != 0) {
            fprintf(debug_file,"%ld -> %ld (w %ld)\n",quant[ox+u][oy+v],
                    dct[ox+u][oy+v],w[u][v]);
          }
          if (dct[ox+u][oy+v]>max) max = dct[ox+u][oy+v];
          if (dct[ox+u][oy+v]<min) min = dct[ox+u][oy+v];
        }
    }

  if (debug_file != 0) {
    fprintf(debug_file,"requantised coefficients are in range [%ld,%ld]\n",
            min,max);
  }
  
  return;
  
}


/*--------------------------------------------------------------------------*/

void init_category_table(long* table) {
  /* initialise lookup table for categories: array index is the number to encode,
     value of the table entry the corresponding category */
  long i,cat;
  long limits[12]={0,1,3,7,15,31,63,127,255,511,1023,2047};
  cat=0;
  for (i=0;i<=2047;i++) {
    if (i>limits[cat]) cat++;
    table[i]=cat;
  }
}

/*--------------------------------------------------------------------------*/

void init_zigzag_table(long* zigzag_x, long* zigzag_y) {
  /* create lookup tables for zig-zag pattern: array index is the number of the
  coefficient in the linear sequence of the zig-zag pattern, starting with 0,
  value is the corresponding coordinate (for x and y tables respectively). */
  long i,j,c;
  int flag;
  c=0;
  i=0;
  j=0;
  zigzag_x[0]=0;
  zigzag_y[0]=0;
  while (i!=7 || j!=7) {
    /* supplement your code here */
    /* HINT: set i and j to the next coordinate in the zig zag pattern.
       You can use case distinctions to do that. */ 
  	flag=0;
  	if((i==0)||(i==7))
  		{
  			j++;
  			c++;
    zigzag_x[c]=i;
    zigzag_y[c]=j;
  //  printf(" \n %ld %ld %ld ",zigzag_x[c],zigzag_y[c],c);
    flag=1;
    if(i==7 && j==7)
    	{
    		/*zigzag_x[c]=i;
    zigzag_y[c]=j;
    printf(" \n %ld %ld %ld ",zigzag_x[c],zigzag_y[c],c);*/
    		break;
    	}
	}
  	if(((j==0)||(j==7))&&(flag!=1))
  		{
  			i++;
  			c++;
    zigzag_x[c]=i;
    zigzag_y[c]=j;
 //   printf(" \n %ld %ld %ld ",zigzag_x[c],zigzag_y[c],c);
}
  	if((i+j)%2!=0)
  	{
  		i++;
  		j--;
  		c++;
    zigzag_x[c]=i;
    zigzag_y[c]=j;
  //  printf(" \n %ld %ld %ld ",zigzag_x[c],zigzag_y[c],c);
  	}
  	if((i+j)%2==0)
  	{
  		i--;
  		j++;
  		c++;
    zigzag_x[c]=i;
    zigzag_y[c]=j;
 //  printf(" \n %ld %ld %ld ",zigzag_x[c],zigzag_y[c],c);
  	}

    
  }
  return;
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
  if (output_file !=0) {
    for (i=0;i<n;i++) {
      set_bit(output_file,c%2);
      c/=2;
    }
  }
  return;
}


/*--------------------------------------------------------------------------*/
/* calulates block DCT of input image/channel */
void block_encode(long  **quant,      /* input quantised DCT coefficients */
                  long nx, long ny,   /* image dimensions */
                  FILE* debug_file,   /* 0 - no output, 
                                         otherwise debug output to file */
                  BFILE *binary_file) {/* file for binary output */

  long u,v,k,l,i;      /* loop variables */
  long N=8;            /* block size */
  long blocks_x;       /* number of blocks in each direction */
  long blocks_y;
  long ox,oy;          /* block offsets */
  long category_lookup[2048]; /* lookup table for categories */
  long zigzag_x[64];   /* x-index for zig-zag traversal of blocks */ 
  long zigzag_y[64];   /* y-index for zig-zag traversal of blocks */
  long last_dc = 0;    /* previously encoded dc coefficient */
  long cat;            /* category */
  long c;              /* number to encode in each category */
  long runlength;      /* run length */
  long symbols;        /* number of symbols to encode */
  long *cache_sym;     /* temporary storage for symbols */
  long *cache_c;       /* temporary storage for numbers */
  long pred_error;     /* prediction error for DC coefficients */

  /* allocate memory */
  alloc_long_vector(&cache_sym,nx*ny);
  alloc_long_vector(&cache_c,nx*ny);
  
  /* determine number of blocks */
  blocks_x = nx/N;
  if ((nx % N) > 0) blocks_x++;
  blocks_y = ny/N;
  if ((ny % N) > 0) blocks_y++;


  if (debug_file != 0) {
    fprintf(debug_file,"ENCODING\n");
    fprintf(debug_file,"========\n");
  }

  /* initialise lookup tables */
  init_category_table((long*)category_lookup);


  init_zigzag_table((long*)zigzag_x,(long*)zigzag_y);
  //printf("zigzagged\n");

  /* initialise symbol counter */
  symbols = 0;
  
  /* Iterate over all blocks and transform them into sequence of symbols */
  for (k=0;k<blocks_x;k++)
    for (l=0;l<blocks_y;l++) {
      ox = k*N+1; oy=l*N+1; /* define block offsets */

      /* print block for debugging */
      if (debug_file != 0) {
        fprintf(debug_file,"Block %ld %ld\n",k,l);
        fprintf(debug_file,"quantised DCT:\n");
        for (v=0; v<N; v++) {
          for (u=0; u<N; u++) {
            fprintf(debug_file,"%ld ",quant[ox+u][oy+v]);
          }
          fprintf(debug_file,"\n");
        }
        fprintf(debug_file,"encoded sequence:\n");
      }

      /* encode DC coefficient of current block */
      pred_error = quant[ox][oy]-last_dc;
      cat = category_lookup[abs(pred_error)];
      /* supplement your code here */
            /* HINT: the array category_lookup has been initialised with
               init_category_table to assist you. Find out what it does
               and use it here. */
      if (pred_error > 0) {
        c = pred_error;
      } else {
        c = (long)pow(2,cat)-1+pred_error;
      }

      /* store DC representation into cache */
      cache_sym[symbols]=cat;
      cache_c[symbols]=c;
      symbols++;
      
      /* write encoded DC coefficient to debug file */
      if (debug_file != 0) {
        fprintf(debug_file,"DC: %ld (diff %ld, last %ld, ",
                cat,pred_error,last_dc);
      }
        write_long_bitwise(c,cat,debug_file,0);
      if (debug_file != 0) {
        fprintf(debug_file,") ");
      }

      last_dc = quant[ox][oy];
      
      /* store AC coefficients */
      /* symbols for AC: symbol:=12*runlength+cat covers 0,...,191 */
      /* EOB: 193, ZRL: 192, both are available as defines */
      /* for symbols without associated c, set cache_c to -1 */
      runlength=0;
      for (i=1;i<64;i++) {
        u=zigzag_x[i]; v=zigzag_y[i];
        if (quant[ox+u][oy+v]==0) {
          runlength++;
        } else {
          cat = category_lookup[abs(quant[ox+u][oy+v])];
          if (quant[ox+u][oy+v] > 0) {
            c = quant[ox+u][oy+v];
          } else {
            c = (long)pow(2,cat)-1+quant[ox+u][oy+v];
          }

          /* handle run lengths > 15 */
          while (runlength > 15) {
            cache_sym[symbols]=ZRL;
            cache_c[symbols]=-1;
            symbols++;
            runlength=0;/* supplement your code here */
            if (debug_file != 0) {
            fprintf(debug_file,"ZRL ");
            }
          }
          
          /* store AC representation into cache */
          cache_sym[symbols]=cat+12*runlength;
          cache_c[symbols]=c;
          symbols++;

          /* write encoded AC coefficient to debug file */
          if (debug_file != 0) {
            fprintf(debug_file,"%ld/%ld ~ %ld (",
                    runlength,cat,cache_sym[symbols-1]);
          }
            write_long_bitwise(c,cat,debug_file,0);
          if (debug_file != 0) {
            fprintf(debug_file,") ");
          }
          runlength=0;
        }
      }

      /* handle end of block (EOB) */
      if (runlength > 0) {
        cache_sym[symbols]=EOB;
        cache_c[symbols]=-1;
        symbols++;
        if (debug_file != 0) {
          fprintf(debug_file,"EOB");
        }
      }

      if (debug_file != 0) {
        fprintf(debug_file,"\n");
      }
    }

  /* encode and store symbols with adaptive arithmetic coding */
  encode_adaptive_wnc(cache_sym,symbols,194,0.3,(long)pow(2,8),0,binary_file);

  /* append category offsets at end of file */
  for (i=0;i<symbols;i++) {
    if (cache_c[i] != -1) {
      /* category is given as symbol mod 12 */
      write_long_bitwise(cache_c[i],cache_sym[i]%12,debug_file,0);
    }
  }

  /* free memory */
  disalloc_long_vector(cache_sym,nx*ny);
  disalloc_long_vector(cache_c,nx*ny);

  return;
  
}

/*--------------------------------------------------------------------------*/
float mse(long ***u, long ***f, long nx, long ny, long nc) {
  long i,j,c;
  long sum, diff;

  sum=0;
  for (c=0;c<nc;c++) 
    for (i=1;i<=nx;i++)
      for (j=1;j<=ny;j++) {
        diff = u[c][i][j]-f[c][i][j];
        sum += diff*diff;
      }
  
  return (double)sum/(double)(nx*ny*(nc+1));
}

/*--------------------------------------------------------------------------*/
void subsample(long **f, /* input: fine resolution */
               long **g, /* output: coarse resolution */
               long nx, long ny, /* fine resolution */
               long factor) {    /* downsampling factor */
  /* subsample by factor in each dimension with simple averaging */
  
  long nx_coarse, ny_coarse, sum, counter, k, l, u, v, ox, oy;
//  int i,j;
  /* determine coarse resolution */
  nx_coarse = nx/factor;
  ny_coarse = ny/factor;
  if ((nx % factor) > 0) nx_coarse++;
  if ((ny % factor) > 0) ny_coarse++;
  
   /* Iterate over all fine resolution blocks corresponding to a coarse
      resolution pixel and average */
  for (k=0;k<nx_coarse;k++)
    for (l=0;l<ny_coarse;l++) {
      ox = k*factor+1; oy=l*factor+1; 
      sum=0; counter=0;
      for (u=0; u<factor; u++)
        for (v=0; v<factor; v++) {
          if ((u+ox <= nx) && (v+oy <=ny)) {
            sum+=f[u+ox][v+oy];
            counter++;
          }
        }
      g[k+1][l+1]=(long)((float)sum/(float)(counter));
    }
  /*  for(i=0;i<16;i++)
      	{
      		printf("\n");
      		for(j=0;j<16;j++)
      		printf(" %ld",f[i][j]);
      	}
      	printf("\n");
      	for(i=0;i<8;i++)
      	{
      		printf("\n");
      		for(j=0;j<8;j++)
      		printf(" %ld",g[i][j]);
      	}
      	printf("\n");*/
}

/*--------------------------------------------------------------------------*/
void upsample(long **f, /* input: coarse resolution */
              long **g, /* output: fine resolution */
              long nx, long ny, /* fine resolution */
              long factor) {    /* upsampling factor */
  /* upsample by factor in each dimension with nearest neighbour */
  
  long nx_coarse, ny_coarse, k, l, u, v, ox, oy;
  //int i,j;
  /* determine coarse resolution */
  nx_coarse = nx/factor;
  ny_coarse = ny/factor;
  if ((nx % factor) > 0) nx_coarse++;
  if ((ny % factor) > 0) ny_coarse++;
  
  /* Iterate over all fine resolution blocks corresponding to a coarse
     resolution pixel and average */
  for (k=0;k<nx_coarse;k++)
    for (l=0;l<ny_coarse;l++) {
      ox = k*factor+1; oy=l*factor+1;
      for (u=0; u<factor; u++)
        for (v=0; v<factor; v++) {
          if (u+ox <= nx && v+oy <=ny) {
            g[u+ox][v+oy]=f[k+1][l+1];
          }
        }

    }
    /*for(i=0;i<8;i++)
      	{
      		printf("\n");
      		for(j=0;j<8;j++)
      		printf(" %ld",f[i][j]);
      	}
      	printf("\n");
      	for(i=0;i<16;i++)
      	{
      		printf("\n");
      		for(j=0;j<16;j++)
      		printf(" %ld",g[i][j]);
      	}
      	printf("\n");*/
    
      	
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
  char tmp_file[1000];         /* string for intermediate filenames */
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
  BFILE* binary_file=0;       /* binary file for compressed bitstring */
  char*  debug_file=0;        /* filename for writing debug information */
  FILE*  dfile=0;             /* file for writing debug information */
  long   q=0;                 /* quantisation parameter */
  long   s=0;                 /* chroma subsampling factor */
  long **tmp_img;             /* temporary image */
  
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
    alloc_long_matrix(&tmp_img,image.nx_ext[0]+2,image.ny_ext[0]+2);

    /* convert to YCbCr space or copy over grey value image */
    if (nc > 1) {
      RGB_to_YCbCr(image.orig_rgb,image.orig_ycbcr,nx[0],ny[0]);
    } else {
      copy_matrix_long(image.orig_rgb[0],image.orig_ycbcr[0],nx[0],ny[0]);
    }
    
    /* perform chroma subsampling */
    if (nc > 1) {
      nx[1]=nx[2]=nx[0]/s;
      ny[1]=ny[2]=ny[0]/s;    
      if ((nx[0] % s) > 0) {nx[1]++;nx[2]++;}
      if ((ny[0] % s) > 0) {ny[1]++;ny[2]++;}
      if (s > 1) {
     /* 	for(i=0;i<nx[1];i++)
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
        subsample(image.orig_ycbcr[1],image.rec_quant[0],nx[0],ny[0],s);
        copy_matrix_long(image.rec_quant[0],image.orig_ycbcr[1],image.nx_ext[1],
                         image.ny_ext[1]);
        subsample(image.orig_ycbcr[2],image.rec_quant[0],nx[0],ny[0],s);
        copy_matrix_long(image.rec_quant[0],image.orig_ycbcr[2],image.nx_ext[1],
                         image.ny_ext[1]);
      }
      printf("Chroma subsampling by factor %ld (%ld x %ld -> %ld x %ld)\n",
             s,nx[0],ny[0],nx[1],ny[1]);
    }
    
    /* extend image dimensions to multiples of block_size */
    for (i=0; i<nc; i++) {
      extend_image(image.orig_ycbcr[i],nx[i],ny[i],image.block_size,
                   image.orig_ycbcr[i]);
    }

    /* open debug file if debug mode is active */
    if (debug_file != 0) {
      dfile = fopen(debug_file,"w");
    }

    /* prepare files for encoding */
    printf("Computing DCT and quantising coefficients\n");
    sprintf(tmp_file,"%s.wnc",output_file);
    binary_file = bfopen(tmp_file,"w");

    /* apply block DCT and encode */
    for (i=0; i<nc; i++) {
      block_DCT(image.orig_ycbcr[i],image.nx_ext[i],image.ny_ext[i],
                image.block_size,image.dct[i]);
      block_quantise(image.dct[i],image.nx_ext[i],image.ny_ext[i],0,
                     image.dct_quant[i]);

      block_encode(image.dct_quant[i],image.nx_ext[i],image.ny_ext[i],
                   dfile,binary_file);          
    }

    /* close binary file */
    bfclose(binary_file);

    /* output image information */
    printf("Resulting compression ratio: %f:1\n\n", 
           get_compression_ratio(input_file,tmp_file));
    
    /* write image data */
    write_comment_string(&image,0,comments);
    for (i=0; i<nc; i++) {
      sprintf(tmp_file,"%s_dct_channel%ld.pgm",output_file,i);
      abs_img(image.dct_quant[i],nx[i],ny[i],tmp_img);
      normalise_to_8bit(tmp_img,nx[i],ny[i],tmp_img);
      write_pgm(tmp_img, nx[i], ny[i],tmp_file, comments);
    }

    /* reconstruct */
    printf("Requantising and compute inverse DCT\n");
    for (i=0; i<nc; i++) {
      block_requantise(image.dct_quant[i],image.nx_ext[i],image.ny_ext[i],0,
                     image.dct_quant[i]);
      block_IDCT(image.dct_quant[i],image.nx_ext[i],image.ny_ext[i],
                 image.block_size,image.rec[i]);
      convert_matrix_int(image.rec[i],image.rec_quant[i],
                        image.nx_ext[i],image.ny_ext[i]);
    }

    /* perform upsampling if downsampling was applied before */
    if (s>1 && nc > 1) {
      upsample(image.rec_quant[1],tmp_img,nx[0],ny[0],s);

      copy_matrix_long(tmp_img,image.rec_quant[1],nx[0],ny[0]);
      upsample(image.rec_quant[2],tmp_img,nx[0],ny[0],s);
      copy_matrix_long(tmp_img,image.rec_quant[2],nx[0],ny[0]);

      printf("Chroma upsampling by factor %ld (%ld x %ld -> %ld x %ld)\n",
             s,nx[1],ny[1],nx[0],ny[0]);
    }

    /* convert back from YCbCr to RGB */
    if (nc>1) {
      YCbCr_to_RGB(image.rec_quant,image.rec_quant,image.nx_ext[0],
                   image.ny_ext[0]);
    }

    printf("Resulting MSE: %f\n",mse(image.orig_rgb,image.rec_quant,nx[0],ny[0],
                                     nc));

    /* write reconstruction */
    if (format==FORMAT_PPM) {
      sprintf(tmp_file,"%s_rec.ppm",output_file);
      write_ppm(image.rec_quant, nx[0], ny[0], tmp_file, comments);
    } else {
      sprintf(tmp_file,"%s_rec.pgm",output_file);
      write_pgm(image.rec_quant[0], nx[0], ny[0], tmp_file, comments);
    }
    
    if (debug_file !=0) {
      fclose(dfile);
     }
     
  } else {
    /* DECOMPRESS *************************************************************/

  }
  
  /* ---- free memory  ---- */
  disalloc_long_matrix(tmp_img,image.nx_ext[0]+2,image.ny_ext[0]+2);
  destroy_image(&image);
  free(program_call);

  return(0);
}


/* -------------------- STORAGE AND CODING ---------------------------------*/

typedef struct _BITBUFFER {
  unsigned char *b; /* array for bytewise content */
  long max_size;    /* maximum size of the buffer (determined by allocation) */
  long byte_pos;    /* current byte position */
  long bit_pos;     /* current bit position (inside current byte) */
} BITBUFFER;

/*--------------------------------------------------------------------------*/

void alloc_bitbuffer

  (BITBUFFER *buffer, /* bitbuffer */
   long n)            /* size */

/* allocates memory for a bitbuffer of size n */
{
long i; /* loop variable */

/* allocate array for individual bytes */
buffer->b = (unsigned char *) malloc(n * sizeof(unsigned char));
if (buffer->b == NULL)
{
   printf("alloc_bitbuffer: not enough memory available\n");
   exit(1);
}

/* set max size, positions, and initialise all bytes with 0 */
buffer->max_size = n;
for (i = 0; i < buffer->max_size; i++)
{
   buffer->b[i] = 0;
}
buffer->byte_pos = 0;
buffer->bit_pos = 0;
return;
} /* alloc_bitbuffer */

/*--------------------------------------------------------------------------*/

long bitbuffer_addbit

  (BITBUFFER *buffer,  /* bitbuffer (output) */
   unsigned char bit)  /* bit to be written */

/* Add a single bit to the buffer at the current bit position.
 * Returns 1 on success, 0 if the buffer is full
 */
{
if ((buffer->byte_pos < buffer->max_size) && (buffer->bit_pos < 8))
{
   if (bit)
      buffer->b[buffer->byte_pos] |= (1 << buffer->bit_pos);
   else
      buffer->b[buffer->byte_pos] &= ~(1 << buffer->bit_pos);
   (buffer->bit_pos)++;
   if (buffer->bit_pos > 7)
   {
      buffer->bit_pos = 0;
      buffer->byte_pos++;
   }
   return 1;
}
return 0; /* adding not successful, memory full */
} /* bitbuffer_addbit */

/*--------------------------------------------------------------------------*/

long bitbuffer_getbit

  (BITBUFFER *buffer)  /* bitbuffer (output) */

/* Get a single bit from the buffer at the current bit position and
 * move the position one bit further.
 * Returns the bit on success, -1 on failure.
 */
{
unsigned char bit;

if ((buffer->byte_pos < buffer->max_size) && (buffer->bit_pos < 8))
{
   bit = (buffer->b[buffer->byte_pos] >> buffer->bit_pos) & 1;
   (buffer->bit_pos)++;
   if (buffer->bit_pos > 7)
   {
      buffer->bit_pos = 0;
      buffer->byte_pos++;
   }
   return bit;
}
else
{
   return -1; /* end of buffer already reached, no more bits to read */
}
} /* bitbuffer_getbit */

/*--------------------------------------------------------------------------*/

void bitbuffer_writefile

  (BITBUFFER *buffer, /* buffer to be written */
   FILE *bitfile)     /* output file (should be open for writing) */

/* Write full bitbuffer to file in a single disk operation. */
{
fwrite(buffer->b, sizeof(unsigned char), buffer->byte_pos + 1, bitfile);
} /* bitbuffer_writefile */

/*--------------------------------------------------------------------------*/

void bitbuffer_loadfile

  (BITBUFFER *buffer, /* output buffer*/
   FILE *bitfile)     /* input file (should be open for reading "rb") */

/* Load file content from current position until end into bitbuffer. */
{
long start_pos; /* initial position of file pointer */
long chunksize; /* size of file chunk to read */


start_pos = ftell(bitfile);
fseek(bitfile, 0, SEEK_END);
chunksize = ftell(bitfile)-start_pos;
fseek(bitfile, start_pos, SEEK_SET);

/*printf("allocating buffer of size %ld\n",chunksize);*/
alloc_bitbuffer(buffer,chunksize);

/*printf("start reading at position %ld, write to %ld",
   start_pos,buffer->byte_pos); */

fread(buffer->b, 1, chunksize, bitfile);
/*printf("%ld bytes loaded\n",chunksize);*/

} /* bitbuffer_loadfile */

/*--------------------------------------------------------------------------*/

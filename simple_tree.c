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
  node->cx = cx;
  node->cy = cy;
  node->nx = nx;
  node->ny = ny;
  node->error = error;
}

/*--------------------------------------------------------------------------*/
void init_tree(const REEDTree *tree, long nx, long ny) {

  long i,n,j;  /* loop variables */
  long cx,cy;  /* corners of subimage */
  long nx_new=0, ny_new=0; /* dimensions of subimage */
  REEDNode node;     /* current node */
  REEDNode parent;     /* parent node */

  /* initialise root node with full image */
  set_node(&(tree->nodes[0][0]),0,1,1,nx,ny,-1);

  for( i = 1; i < (long)(tree->level); i++ ) {
    n = pow( 2.0, i );
      for( j = 0; j < n; j++ ) {
        node = tree->nodes[i][j];
        parent = *(node.parent);
        /* split subimage in its largest dimension */
        if ((j%2)==0) { /* left child */
          cx=parent.cx; cy=parent.cy;
          if (parent.nx>=parent.ny) {
            nx_new=parent.nx/2+1;
            ny_new=parent.ny;
          } else {
            nx_new = parent.nx;
            ny_new = parent.ny/2+1;
          }
        } else { /* right child */
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
        set_node(&(tree->nodes[i][j]),-1,cx,cy,nx_new,ny_new,-1);
      }
    }
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
      FILE *fp)                  /* bitfile  to write to */
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
      FILE *fp)              /* bitfile  to write to */
    /* write reed tree to bitfile,
       return 0 on success, 1 else. */
{
    unsigned long i = 0, j = 0;     /* loop variables */
    unsigned long n = 0;            /* number of nodes on current level */
    long bit = 0;                    /* single bit, being loaded from bitfile */
    unsigned char byte;
    long bit_pos;

    /* validate pointer */
    if( tree == NULL || fp == NULL )
       return 1;

    /* 1. header */
    fread (&byte, sizeof(char), 1, fp);
    tree->min_level = (long)byte;
    fread (&byte, sizeof(char), 1, fp);
    tree->max_level = (long)byte;

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


    return 0;
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

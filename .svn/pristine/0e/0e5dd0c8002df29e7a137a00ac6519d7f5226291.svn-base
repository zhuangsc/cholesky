#include "chol_utils.h"


#include <stdlib.h>


#if 1
void print_matrix_blocked( FILE* str, char * name, int m, int mr, int n, int nr, int br, int bc, double * A ) {
	int i;
	for( i = 0; i < m; i++ ) {
		int j;
		for( j = 0; j < n; j++ ) {
			double e = A[ (j/bc)*mr*bc + (j%bc)*br + (i/br)*br*bc + (i%br) ];
			fprintf(str, "%s(%d,%d) = %.14f;\n", name, i+1, j+1,e);
			//printf("%s(%d,%d) = %.14f;\n", name, i+1, j+1,e);
			//if(e==0) printf("%s: zero at %i,%i\n",name,i,j);
		}
	}
}



double * canonical_from_hierarchical ( char * name, int m, int mr, int n, int nr, int br, int bc, double * A ) {
	double * Atmp = (double *) malloc( m * n * sizeof( double ) );
	if (Atmp == NULL) { 
		perror("Cannot allocate memory"); 
		exit(1); 
	}

	int i;
	for( i = 0; i < m; i++ ) {
		int j;
		for( j = 0; j < n; j++ ) {
			int offset = (j/bc)*mr*bc + (j%bc)*br + (i/br)*br*bc + (i%br);
			Atmp[j*m+i] = A[ offset ];
			//printf("checking (%i,%i) offs %i : %.2f\n", i, j, offset, A[offset]);
			//double *ap= &A[ (j/bc)*mr*bc + (j%bc)*br + (i/br)*br*bc + (i%br) ] ;
			//if(Atmp[j*m+i]==0) printf("found zero at get of %s\n",name);
    		}
  	}
  
#if 0
  for( i = m; i < mr; i++ ) 
  {
    int j;
    for( j = n; j < nr; j++ ) 
    {
      double padd = A[ (j/ts)*mr*ts + (j%ts)*ts + (i/ts)*ts*ts + (i%ts) ] ;
      if(padd!=0.0)
      {
	fprintf(stderr,"bad pad at (%i,%i) for %s\n",i,j,name);
      }
    }
  }
#endif

	return Atmp;
}


void print_matrix_canonical( char *name, int m, int n, double *A, int lda ) 
{
  int i;
  for (i=0; i<m; i++)
  {
    int j;
    for (j=0; j<n; j++)
    {
      printf("%s(%d,%d) = %22.16g;\n", name, i+1,j+1, A[j*lda+i] );
    }
  }
}
#endif

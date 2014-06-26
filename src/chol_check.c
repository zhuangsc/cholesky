#include "chol_check.h"

#include <math.h>
#include <stdlib.h>

#include "chol_utils.h"
#include "blas.h"




#if 0
void qrca_zerocheck(int m, int mr, int n, int nr, int ts, double *A, char *file, int line)
{
  double *R = get_canonical_from_hierarchical( "A", m, mr, n, nr, ts, A);

  int j;
  for(j=0;j<n;j++)
  {
    int i;
    for(i=0;i<m;i++)
    {
      if(R[j*m+i]==0.0)
	printf("zero at %s:%i\n",file,line);
    }
  }

  free(R);
}
#endif

// Perform equivalent to Octave's command:
// norm(triu(abs(M1)-abs(M2)),1)/norm(M1,1);
double __check ( double * M1, double * M2, int m) {
	double *diffmat  = (double *) calloc( (m*m), sizeof(double) );
	if (diffmat == NULL) { 
		perror("Error Allocating diffmat"); 
		exit(1); 
	}
	double *sumvec   = (double *) calloc( (m),   sizeof(double) );
  	if (sumvec == NULL) { 
		perror("Error Allocating sumvec"); 
		exit(1); 
	}
	double *sumvecM1 = (double *) calloc( (m),   sizeof(double) );
  	if (sumvecM1 == NULL) { 
		perror("Error Allocating sumvecM1"); 
		exit(1); 
	}


  	int j;
	for (j=0; j<m; j++) {
		int i;
		for (i=j; i<m; i++) {
      			double tmp = fabs(M1[j*m+i]);
      			sumvecM1[j] += tmp;
      			tmp -= fabs(M2[j*m+i]);
			//printf("dif (%i,%i) %x\n", i, j, &diffmat[j*m+i]);
      			diffmat[j*m+i] = tmp;
	      		sumvec[j] += fabs(tmp);
    		}
  	}

	double norm1=0.0, norm1M1=0.0;
	for (j=0; j<m; j++) {
		norm1= max(norm1,sumvec[j]);
		norm1M1 = max(norm1M1,sumvecM1[j]);
	}
	norm1 /= norm1M1;

	free(sumvec);
	free(sumvecM1);
	free(diffmat);

	return(norm1);
}


int chol_check(int check, int m, int mr, int ts, double **Ah, double *Rcheck) {
	if( check > 0 ) {
#if USE_PRL
    		double *R = canonical_from_hierarchical( "R", m, mr, m, mr, mr, ts, Ah[ 0 ] );
#else
    		double *R = canonical_from_hierarchical( "R", m, mr, m, mr, ts, ts, Ah[ 0 ] );
#endif
		printf( "calling dpotrf for residual checking...");
		fflush(0);

		int IINFO;
		dpotrf_("Lower", &m, Rcheck, &m, &IINFO);
		printf( "...done (IINFO: %d)\n", IINFO);
		if(IINFO < 0) {
			perror("Error in dpotrf\n");
		}

		double norm_diff = __check( Rcheck, R, m);
		printf( "norm of difference: %22.16g\n", norm_diff);

		free(Rcheck);

		int ret=(norm_diff > 0.0001)?1:0;

		return ret;
	}

	return 0;
}

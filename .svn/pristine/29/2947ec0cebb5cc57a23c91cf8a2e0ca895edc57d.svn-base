#include "genmat.h"

#include "fptype.h"

#include <stdlib.h>
#include <stdio.h>


// generates an SPD matrix in blocked layout in place
void GENMAT(int n, int tn, int nleft, int b, fp_t *A) {
	fp_t *Ain = &A[0];

	int j;
	for (j = 0; j < tn; j++ ) {
		int jend=j==tn-1;
		int i;
		for( i = 0; i < tn; i++ ) {
			int iend=i==tn-1;
			int jj; 
			for( jj=0 ; jj < b; jj++ ) {
				int ii;
				for ( ii=0 ; ii < b; ii++ ) {
  					int skipright = iend?ii>=b-nleft:0;
  					int skipbotn = jend?jj>=b-nleft:0;
	  
  					if ( !skipright && !skipbotn ) {
    						fp_t dran = drand48();
    						if ( dran == 0.0 )  printf("generated 0\n");
							if( i==j && ii ==jj ) {
								dran += n;
							}
    						*Ain = dran;

    						//printf("A(%i,%i)=%.4f\n",i*b+ii+1,j*b+jj+1,dran);
  					} else {
						if ( i==j & ii==jj ) {
							*Ain=1.00;
						}
					}
					Ain++;
				}
			}
		}
  	}
}


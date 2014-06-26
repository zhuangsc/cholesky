#include "chol_config.h"

#include "hb.h"
#include "iohb.h"
#include "hbpad.h"

#include <sys/time.h>


extern int check;
extern int m;
extern int morig;
extern int ts;
extern int bs;
extern int mr;
extern int reps;
extern void *Ahb;
extern int format;

int chol_config(int argc, char *argv[]) {
	if ( argc  < 4 ) {
		fprintf(stderr, "usage %s HBfile b CSC(0)|CSR(1) [rep]\n", argv[0]);
		return 1;
	}

	bs = ts = atoi(argv[2]);
	format = atoi(argv[3]);
	reps = 1;
	if ( argc > 4 ) {
		reps = atoi(argv[4]);
	}

	hbmat_t *lAhb = Ahb = malloc(sizeof(hbmat_t));

	struct timeval start, stop;
	gettimeofday(&start, NULL);
	readHB_newmat_double(argv[1], &(lAhb->m), &(lAhb->n), &(lAhb->elemc), &(lAhb->vptr), &(lAhb->vpos), (double **)&(lAhb->vval));
	gettimeofday(&stop, NULL);
	double elapsed = (stop.tv_sec - start.tv_sec) * 1e6 + stop.tv_usec - start.tv_usec;
	printf("info: read %s in %f ms\n", argv[1], elapsed / 1e3);

	mr = m = lAhb->m;
	if ( m % bs != 0 ) {
		mr = (( m + bs - 1 ) / bs ) * bs;
		Ahb = csc_pad(Ahb, mr, OFFS_F);
		if ( Ahb == NULL ) {
			return 1;
		}
	}

	check = 0;

	return 0;
}

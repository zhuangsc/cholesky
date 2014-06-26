#include "chol_config.h"

#include <stdlib.h>
#include <stdio.h>


extern int m;
extern int ts;
extern int bs;
extern int reps;
extern int check;
extern int mt;
extern int mr;
extern int mtleft;


int chol_config(int argc, char *argv[]) {
	if (argc < 3 ) {
		printf("usage: %s n ts [bs] [reps] [check]\n",argv[0]);
		return 1;
	}

	m = atoi(argv[1]);
	ts = atoi(argv[2]);

	bs = ts;
	reps = 1;
	check = 0;

	if ( argc > 3 ) {
		bs = atoi(argv[3]);
		if ( argc > 4 ) {
			reps = atoi(argv[4]);
			if ( argc > 5 ) {
  				check = atoi(argv[5]);
			}
		}
  	}

	mt = (m+ts-1) / ts;
  	mr = mt * ts;
  	mtleft = mr - m; 

  	printf("A blocks %ix%i elems %ix%i->%ix%i\n",mt,mt,m,m,mr,mr);

	return 0;
}


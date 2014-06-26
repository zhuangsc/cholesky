#ifndef __VECTOR_H__
#define __VECTOR_H__


#include <stddef.h>
#include <malloc.h>
#include <stdlib.h>


#define VECT_DEFSIZE	128

typedef union struc_velem {
	void *p;
	int i;
	double d;
} vel_t;

typedef struct struct_vect {
	vel_t* elem;
	size_t elemc;
	size_t size;
} vector_t;


// debug
static inline void __attribute__((always_inline)) vector_printi(vector_t *v) {
	vel_t *elem = v->elem;
	size_t elemc = v->elemc;

	int i;
	for ( i=0; i<elemc; i++ ) {
		printf("%i ", elem[i].i);
	}	

	if ( elemc ) {
		printf("\n");
	}
}

static inline void __attribute__((always_inline)) vector_printd(vector_t *v) {
	vel_t *elem = v->elem;
	size_t elemc = v->elemc;

	int i;
	for ( i=0; i<elemc; i++ ) {
		printf("%.4lf ", elem[i].d);
	}	

	if ( elemc ) {
		printf("\n");
	}
}


static inline int __attribute__((always_inline)) vector_contains(vector_t *v, vel_t obj) {
	size_t elemc = v->elemc;
	vel_t *elem = v->elem;
	
	int i;
	for ( i=0; i<elemc; i++ ) {
		if ( elem[i].i == obj.i ) {
			return 1;
		}
	}
	
	return 0;
}


// creation / destruction 
static inline vector_t* __attribute__((always_inline)) vector_create() {
	vector_t *v = (vector_t*) malloc(sizeof(vector_t));
	vel_t *elem = (vel_t*) malloc(sizeof(vel_t) * VECT_DEFSIZE);

	v->elem = elem;
	v->elemc = 0;
	v->size = VECT_DEFSIZE;

	return v;
}

static inline vector_t* __attribute__((always_inline)) vector_create_size(int in_size) {
	vector_t *v = (vector_t*) malloc(sizeof(vector_t));
	vel_t *elem = (vel_t*) malloc(sizeof(vel_t) * in_size);

	v->elem = elem;
	v->elemc = 0;
	v->size = in_size;

	return v;
}


static inline void __attribute__((always_inline)) vector_free(vector_t *v) {
	if ( v->elem != NULL ) {
		free(v->elem);
	}
	free(v);
}


static inline vel_t __attribute__((always_inline)) vector_get(vector_t *v, size_t i) {
	return v->elem[i];
}


static inline void __attribute__((always_inline)) vector_insertat(vector_t *v, vel_t obj, size_t i) {
	if ( i >= v->size ) {
		size_t newsize = v->size + VECT_DEFSIZE;
		v->elem = realloc(v->elem, sizeof(vel_t) * newsize);
		v->size = newsize;
	}

	size_t s = i+1;
	if ( s > v->elemc ) {
		v->elemc = s;
	}

	v->elem[i] = obj;
}

static inline void __attribute__((always_inline)) vector_insert(vector_t *v, vel_t obj) {
	vector_insertat(v, obj, v->elemc);
}


static inline int __attribute__((always_inline)) vector_insert_t(vector_t *v, vel_t obj) {
	if ( ! vector_contains(v, obj) ) {
		vector_insertat(v, obj, v->elemc);
		return 1;
	}

	return 0;
}

static inline int __attribute__((always_inline)) vector_contains_partial(vector_t *v, vel_t obj, int from) {
	size_t elemc = v->elemc;
	vel_t *elem = v->elem;
	
	int i;
	for ( i=from; i<elemc; i++ ) {
		if ( elem[i].i == obj.i ) {
			return 1;
		}
	}
	
	return 0;
}

static inline int __attribute__((always_inline)) vector_insert_t_partial(vector_t *v, vel_t obj, int from) {
	if ( ! vector_contains_partial(v, obj, from) ) {
		vector_insertat(v, obj, v->elemc);
		return 1;
	}

	return 0;
}

static inline vector_t* __attribute__((always_inline)) vector_appendi(vector_t *v, int *objs, size_t size) {
	int i;
	for ( i=0; i< size; i++ ) {
		vel_t e;
		e.i = objs[i];
		vector_insert(v, e);
	}

	return v;
}


static inline vector_t* __attribute__((always_inline)) vector_append(vector_t *v, vel_t *objs, size_t size) {
	int i;
	for ( i=0; i< size; i++ ) {
		vector_insert(v, objs[i]);
	}

	return v;
}


static inline void __attribute__((always_inline)) vector_clear(vector_t *v) {
	v->elemc = 0;
}


static inline void __attribute__((always_inline)) vector_rm(vector_t *v, vel_t obj) {
	vel_t* elem = v->elem;
	size_t elemc = v->elemc;

	int i=-1;
	int found = 0;
	while ( !found ) {
		++i;
		found = elem[i].i==obj.i;
	}

	int j;
	for ( j=i; j<elemc; j++ ) {
		elem[j] = elem[j+1];
	}

	if ( i == elemc - 1) {
		--v->elemc;
	}
}

// sorting
static inline int vector_compare_int(const void *a, const void *b){
	const vel_t *da = (const vel_t *)a; 
	const vel_t *db = (const vel_t *)b; 
	const vel_t daa = *da;
	const vel_t dbb = *db;

	return (daa.i > dbb.i) - (daa.i < dbb.i);
}

static inline void __attribute__((always_inline)) vector_qsorti(vector_t *v){
	vel_t *elem = v->elem;
	size_t elemc = v->elemc;
	qsort(elem, elemc, sizeof(vel_t), vector_compare_int);
}

static inline void __attribute__((always_inline)) vector_qsorti_partial(vector_t *v, int from){
	vel_t *elem = v->elem;
	size_t elemc = v->elemc;
	qsort(&(elem[from]), elemc-from, sizeof(vel_t), vector_compare_int);
}

// conversion
static inline int* __attribute__((always_inline)) vector2int(vector_t *v) {
	vel_t *elem = v->elem;
	size_t elemc = v->elemc;
	int *intelem = (int*) malloc( v->elemc * sizeof(int) );

	int i;
	for ( i=0; i < elemc ; i++ ) {
		intelem[i] = elem[i].i;
	}
	
	vector_free(v);

	return intelem;
}

static inline double* __attribute__((always_inline)) vector2double(vector_t *v) {
	vel_t *elem = v->elem;
	size_t elemc = v->elemc;
	double *delem = (double*) malloc( v->elemc * sizeof(double) );

	int i;
	for ( i=0; i < elemc ; i++ ) {
		delem[i] = elem[i].d;
	}
	
	vector_free(v);

	return delem;
}

static inline int* __attribute__((always_inline)) vector2int_nf(vector_t *v) {
	vel_t *elem = v->elem;
	size_t elemc = v->elemc;
	int *intelem = (int*) malloc( v->elemc * sizeof(int) );

	int i;
	for ( i=0; i < elemc ; i++ ) {
		intelem[i] = elem[i].i;
	}
	
//	vector_free(v);

	return intelem;
}

static inline double* __attribute__((always_inline)) vector2double_nf(vector_t *v) {
	vel_t *elem = v->elem;
	size_t elemc = v->elemc;
	double *delem = (double*) malloc( v->elemc * sizeof(double) );

	int i;
	for ( i=0; i < elemc ; i++ ) {
		delem[i] = elem[i].d;
	}
	
//	vector_free(v);

	return delem;
}



//Geschwindigkeit

typedef struct struct_vect_int {
	int* elem;
	size_t elemc;
	size_t size;
} vector_int;

static inline void __attribute__((always_inline)) vector_int_clear(vector_int *v) {
	v->elemc = 0;
}

static inline vector_int* __attribute__((always_inline)) vector_int_create(int *array, int vect_size) {
	vector_int *v = malloc(sizeof(vector_int));
	int *elem = array;

	v->elem = elem;
	v->elemc = 0;
	v->size = vect_size;
//	v->size = VECT_DEFSIZE;

	return v;
}

static inline void __attribute__((always_inline)) vector_int_insertat(vector_int *v, int obj, size_t i) {
	if ( i >= v->size ) {
		size_t newsize = v->size + VECT_DEFSIZE;
		v->elem = realloc(v->elem, sizeof(int) * newsize);
		v->size = newsize;
	}

	size_t s = i+1;
	if ( s > v->elemc ) {
		v->elemc = s;
	}

	v->elem[i] = obj;
}

static inline void __attribute__((always_inline)) vector_int_insert(vector_int *v, int obj) {
	vector_int_insertat(v, obj, v->elemc);
}

static inline int __attribute__((always_inline)) vector_int_contains_partial(vector_int *v, int obj, int from) {
	size_t elemc = v->elemc;
	int *elem = v->elem;
	
	int i;
	for ( i=from; i<elemc; i++ ) {
		if ( elem[i] == obj ) {
			return 1;
		}
	}
	
	return 0;
}

static inline int __attribute__((always_inline)) vector_int_insert_t_partial(vector_int *v, int obj, int from) {
	if ( ! vector_int_contains_partial(v, obj, from) ) {
		vector_int_insertat(v, obj, v->elemc);
		return 1;
	}

	return 0;
}

static inline int vector_int_compare(const void *a, const void *b){
	const int *da = (const int *)a; 
	const int *db = (const int *)b; 
	const int daa = *da;
	const int dbb = *db;

	return (daa > dbb) - (daa < dbb);
}

static inline void __attribute__((always_inline)) vector_int_qsorti_partial(vector_int *v, int from){
	int *elem = v->elem;
	size_t elemc = v->elemc;
	qsort(&(elem[from]), elemc-from, sizeof(int), vector_int_compare);
}

static inline void __attribute__((always_inline)) vector_int_print(vector_int *v) {
	int *elem = v->elem;
	size_t elemc = v->elemc;

	int i;
	for ( i=0; i<elemc; i++ ) {
		printf("%i ", elem[i]);
	}	

	if ( elemc ) {
		printf("\n");
	}
}




typedef struct struct_vect_double {
	double* elem;
	size_t elemc;
	size_t size;
} vector_double;

static inline void __attribute__((always_inline)) vector_double_clear(vector_double *v) {
	v->elemc = 0;
}

static inline vector_double* __attribute__((always_inline)) vector_double_create(double *array, int vect_size) {
	vector_double *v = malloc(sizeof(vector_double));
	double *elem = array;

	v->elem = elem;
	v->elemc = 0;
	v->size = vect_size;

	return v;
}

static inline void __attribute__((always_inline)) vector_double_insertat(vector_double *v, double obj, size_t i) {
	if ( i >= v->size ) {
		size_t newsize = v->size + VECT_DEFSIZE;
		v->elem = realloc(v->elem, sizeof(double) * newsize);
		v->size = newsize;
	}

	size_t s = i+1;
	if ( s > v->elemc ) {
		v->elemc = s;
	}

	v->elem[i] = obj;
}

static inline void __attribute__((always_inline)) vector_double_insert(vector_double *v, double obj) {
	vector_double_insertat(v, obj, v->elemc);
}

static inline int __attribute__((always_inline)) vector_double_contains_partial(vector_double *v, double obj, int from) {
	size_t elemc = v->elemc;
	double *elem = v->elem;
	
	int i;
	for ( i=from; i<elemc; i++ ) {
		if ( elem[i] == obj ) {
			return 1;
		}
	}
	
	return 0;
}

static inline int __attribute__((always_inline)) vector_double_insert_t_partial(vector_double *v, double obj, int from) {
	if ( ! vector_double_contains_partial(v, obj, from) ) {
		vector_double_insertat(v, obj, v->elemc);
		return 1;
	}

	return 0;
}

static inline int vector_double_compare(const void *a, const void *b){
	const double *da = (const double *)a; 
	const double *db = (const double *)b; 
	const double daa = *da;
	const double dbb = *db;

	return (daa > dbb) - (daa < dbb);
}

static inline void __attribute__((always_inline)) vector_double_print(vector_double *v) {
	double *elem = v->elem;
	size_t elemc = v->elemc;

	int i;
	for ( i=0; i<elemc; i++ ) {
		printf("%.4lf ", elem[i]);
	}	

	if ( elemc ) {
		printf("\n");
	}
}

#endif // __VECTOR_H__

#ifdef __cplusplus
extern "C" {
#endif


#ifdef LIBOMPSS_BUILDING
#define LIBOMPSS_DLL_EXPORTED __attribute__((__visibility__("default")))
#else
#define LIBOMPSS_DLL_EXPORTED
#endif


#include "hb.h"


/* 	A: (in)		SPD matrix in HB, order n
	b: (in)		block size
	work: (in)	work array, size n
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED hbmat_t* ompss_csr_dchol_ll(int b, hbmat_t *A, int *work);  

extern LIBOMPSS_DLL_EXPORTED hbmat_t* ompss_csc_dchol_ll(int b, hbmat_t *A, int *work);


/* 	A: (in)		SPD matrix in hypermatrix form
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_schol_hll(int n, int b, int t, float **Ah);

/* 	A: (in)		SPD matrix in hypermatrix form
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_hll(int n, int b, int t, double **Ah);

/* 	A: (in)		SPD matrix in hypermatrix form
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_schol_hrl(int n, int b, int t, float **Ah);

/* 	A: (in)		SPD matrix in hypermatrix form
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_hrl(int n, int b, int t, double **Ah);




/* 	A: (in)		SPD matrix 
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_ll(int n, int b, int t, double *A);


/* 	A: (in)		SPD matrix in CSR format
	L: (out)	lower-triangular factor */
//extern LIBOMPSS_DLL_EXPORTED hbmat_t* ompss_csr_dchol_ll(int b, hbmat_t *A);


/* 	A: (in)		SPD matrix in CSR format
	L: (out)	lower-triangular factor */
//extern LIBOMPSS_DLL_EXPORTED hbmat_t* ompss_csc_dchol_ll(int b, hbmat_t *A);



extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_nhll(int mt, int b, int t, double **Ah);

extern LIBOMPSS_DLL_EXPORTED int ompss_schol_nhll(int mt, int b, int t, float **Ah);

extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_nhrl(int mt, int b, int t, double **Ah);

extern LIBOMPSS_DLL_EXPORTED int ompss_schol_nhrl(int mt, int b, int t, float **Ah);


#ifdef __cplusplus
}
#endif

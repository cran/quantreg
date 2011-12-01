/*--- Declarations of routines implemented in
 *  ./sparskit2.f
 *  ------------
 * also called by  the routines in srqfn.c and srqfnc.c
 *				~~~~~~~~~~     ~~~~~~~~
 */

#include <R.h>

/* SparseM -- ./sparskit2.f : */
int F77_NAME(aplb)(int *, int *, int *, double *, int *, int *,
		   double *, int *, int *, double *, int *, int *,
		   int *, int *, int *);
int F77_NAME(amub)(int *, int *, int *, double *,
		   int *, int *, double *, int *, int *,
		   double *, int *, int *, int *, int *, int *);
int F77_NAME(amux)(int *, double *, double *,
		   double *, int *, int *);
int F77_NAME(amudia)(int *, int *, double *,
		     int *, int *, double *, double *, int *, int *);

int F77_NAME(extract)(double *, int *, int *,
		      double *, int *, int *, int *, int *, int *);


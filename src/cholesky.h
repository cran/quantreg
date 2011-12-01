/*--- Declarations of routines implemented in
 *  ./cholesky.c
 *  ------------
 * also called by  chlfct_() in ./chlfct.c
 *				~~~~~~~~~~
 */

#include <R.h>

/* heritage from  f2c : */
typedef int /* Unknown procedure type */ (*U_fp)();
typedef int /* Subroutine */		 (*S_fp)();

int F77_NAME(mmpy1)(int *m, int *n, int *q, int *xpnt,
		    double *x, double *y, int *ldy);
int F77_NAME(mmpy2)(int *m, int *n, int *q, int *xpnt,
		    double *x, double *y, int *ldy);
int F77_NAME(mmpy4)(int *m, int *n, int *q, int *xpnt,
		    double *x, double *y, int *ldy);
int F77_NAME(mmpy8)(int *m, int *n, int *q, int *xpnt,
		    double *x, double *y, int *ldy);

int F77_NAME(inpnv)(int *neqns, int *xadjf, int *adjf,
		    double *anzf, int *perm, int *invp, int *nsuper,
		    int *xsuper, int *xlindx, int *lindx, int *xlnz,
		    double *lnz, int *offset);

int F77_NAME(smxpy1)(int *m, int *n, double *y, int *apnt, double *a);
int F77_NAME(smxpy2)(int *m, int *n, double *y, int *apnt, double *a);
int F77_NAME(smxpy4)(int *m, int *n, double *y, int *apnt, double *a);
int F77_NAME(smxpy8)(int *m, int *n, double *y, int *apnt, double *a);

int F77_NAME(blkfct)(int *neqns, int *nsuper, int *xsuper,
		     int *snode, int *split, int *xlindx, int *lindx,
		     int *xlnz, double *lnz, int *iwsiz, int *iwork,
		     int *tmpsiz, double *tmpvec, int *iflag,
		     U_fp mmpyn, U_fp smxpy);

/* TODO(MM): should call blkfc2() from "above" and drop blkfct() completely */
void
F77_NAME(blkfc2)(int *nsuper, int *xsuper, int *snode, int *split,
		 int *xlindx, int *lindx, int *xlnz, double *lnz,
		 int *link, int *length, int *indmap, int *relind,
		 int *tmpsiz, double *tmpvec, int *iflag,
		 U_fp mmpyn, U_fp smxpy);

int F77_NAME(bfinit)(int *neqns, int *nsuper, int *xsuper,
		     int *snode, int *xlindx, int *lindx, int *cachsz,
		     int *tmpsiz, int *split);

void F77_NAME(ordmmd)(int *neqns, int *xadj, int *adjncy,
		      int *invp, int *perm, int *iwsiz, int *iwork,
		      int *nofsub, int *iflag);

void
F77_NAME(sfinit)(int *neqns, int *nnza, int *xadj,
		 int *adjncy, int *perm, int *invp, int *colcnt,
		 int *nnzl, int *nsub, int *nsuper, int *snode,
		 int *xsuper, int *iwsiz, int *iwork, int *iflag);

int F77_NAME(symfct)(int *neqns, int *adjlen, int *xadj,
		     int *adjncy, int *perm, int *invp, int *colcnt,
		     int *nsuper, int *xsuper, int *snode, int *nofsub,
		     int *xlindx, int *lindx, int *xlnz, int *iwsiz,
		     int *iwork, int *flag);

/* TODO(MM): should call symfc2() from "above" and drop symfct() completely */
void
F77_NAME(symfc2)(int *, int *, int *, int *, int *, int *, int *, int *, int *,
		 int *, int *, int *, int *, int *, int *, int *, int *, int *);

int F77_NAME(blkslv)(int *nsuper, int *xsuper, int *xlindx, int *lindx,
		     int *xlnz, double *lnz, double *rhs);

double F77_NAME(gtimer)(void);

/* ./chlfct.c : */
void F77_NAME(chlfct)(int *m, int *xlindx, int *lindx,
		      int *invp, int *perm, int *iwork, int *nnzdsub,
		      int *jdsub, int *colcnt, int *nsuper, int *snode,
		      int *xsuper, int *nnzlmax, int *nsubmax, int *xlnz,
		      double *lnz, int *id, int *jd, double *d,
		      int *cachsz, int *tmpmax, int *level, double *tmpvec,
		      int *split,
		      int *ierr, /* on error return, sets error code in {3:12},
				    docu., see e.g., ./srqfn.c */
		      int *it, /* it <= 1: do initializations;
				*  otherwise: no init.*/
		      double *timewd);/* of length >= 6:
				       * will return timing information */


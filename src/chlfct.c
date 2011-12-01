/*-*- mode: C; kept-old-versions: 12;  kept-new-versions: 20; -*-
 *
 * chlfct.f -- translated by f2c (version 20031025) and by
 * $Id: f2c-clean,v 1.10 2002/03/28 16:37:27 maechler Exp $
 * plus extended manual code cleaning by Martin Maechler, ETH Zurich
 */

#include <math.h>

/* Subroutine to perform Cholesky factorization
 *			 ======================
 * based on subroutines in ./cholesky.c , declared here :
*/
#include "cholesky.h"

void
F77_SUB(chlfct)(int *m, int *xlindx, int *lindx,
		int *invp, int *perm, int *iwork, int *nnzdsub,
		int *jdsub, int *colcnt, int *nsuper, int *snode,
		int *xsuper, int *nnzlmax, int *nsubmax, int *xlnz,
		double *lnz, int *id, int *jd, double *d,
		int *cachsz, int *tmpmax, int *level, double *tmpvec,
		int *split,
		int *ierr, /* on error return, sets error code in {3:12},
			      docu., see e.g., ./srqfn.c */
		int *it,   /* it <= 1: do initializations;  otherwise: no init.*/
		double *timewd)/* of length >= 6: will return timing information */
/*  as follows: [1] ordmmd()
 *		[2] sfinit() - initialization for sym.fact.
 *		[3] symfct() - symbolic factorization
 *		[4] inpnv()  - Input numerical values into data structures of L
 *		[5] bfinit() - Initialization for block factorization
 *		[6] blkfct() - Numerical factorization
 */
{
    /* Local variables */
    int i, nsub, nnzl, iwsiz, tmpsiz;
    double timbeg;

    --timewd; /* so we can use 1-indexing below */

/* Save the matrix structure from jdsub(m+2:nzzd+1),jdsub(1:m+1)
   to lindx and xlindx because the matrix structure is destroyed by the
   minimum degree ordering routine */

    for (i = 0; i < (*m + 1); ++i)
	xlindx[i] = jdsub[i];
    for (i = 0; i < *nnzdsub; ++i)
	lindx[i] = jdsub[*m + 1 + i];


    iwsiz = *m << 2;

    if (*it <= 1) { /* do the initializations : */

	/* Reorder the matrix using minimum degree ordering routine */

	timbeg = F77_CALL(gtimer)();
	F77_CALL(ordmmd)(m, xlindx, lindx, invp, perm,
			 &iwsiz, iwork, &nsub, ierr);
	timewd[1] += F77_CALL(gtimer)() - timbeg;

	if (*ierr == -1)	{ *ierr = 3; return; }

	/* Call sfinit: Symbolic factorization initialization
	   to compute supernode partition and storage requirements
	   for symbolic factorization. New ordering is a postordering
	   of the nodal elimination tree */

	iwsiz = *m * 7 + 3;
	timbeg = F77_CALL(gtimer)();
	F77_CALL(sfinit)(m, nnzdsub, jdsub, &jdsub[*m + 1], perm, invp,
			 colcnt, &nnzl, &nsub, nsuper, snode, xsuper,
			 &iwsiz, iwork, ierr);
	timewd[2] += F77_CALL(gtimer)() - timbeg;

	if (*ierr == -1)     { *ierr = 4; return; }
	if (nnzl > *nnzlmax) { *ierr = 5; return; }
	if (nsub > *nsubmax) { *ierr = 6; return; }
    }

/* Call symfct: Perform supernodal symbolic factorization */

    /* iwsiz = nsuper + 2 * m + 1; */
    timbeg = F77_CALL(gtimer)();
    F77_CALL(symfct)(m, nnzdsub, jdsub, &jdsub[*m + 1], perm, invp,
		     colcnt, nsuper, xsuper, snode, &nsub,
		     xlindx, lindx, xlnz, &iwsiz, iwork, ierr);
    timewd[3] += F77_CALL(gtimer)() - timbeg;

    if (*ierr == -1)	{ *ierr = 7; return; }
    if (*ierr == -2)	{ *ierr = 8; return; }

/* Call inpnv: Input numerical values into data structures of L */

    timbeg = F77_CALL(gtimer)();
    F77_CALL(inpnv)(m, id, jd, d, perm, invp, nsuper, xsuper,
		    xlindx, lindx, xlnz, lnz, iwork);
    timewd[4] += F77_CALL(gtimer)() - timbeg;

/* Call bfinit: Initialization for block factorization */

    timbeg = F77_CALL(gtimer)();
    F77_CALL(bfinit)(m, nsuper, xsuper, snode, xlindx, lindx,
		     cachsz, &tmpsiz, split);
    timewd[5] += F77_CALL(gtimer)() - timbeg;

    if (tmpsiz > *tmpmax) { *ierr = 9; return; }

/* Call blkfct: Numerical factorization */

    iwsiz = (*m << 1) + (*nsuper << 1);
    timbeg = F77_CALL(gtimer)();
    if (*level == 1) {
	F77_CALL(blkfct)(m, nsuper, xsuper, snode, split,  xlindx, lindx,
			 xlnz, lnz, &iwsiz, iwork, &tmpsiz, tmpvec, ierr,
			 (U_fp) F77_NAME(mmpy1), (U_fp) F77_NAME(smxpy1));
    } else if (*level == 2) {
	F77_CALL(blkfct)(m, nsuper, xsuper, snode, split,  xlindx, lindx,
			 xlnz, lnz, &iwsiz, iwork, &tmpsiz, tmpvec, ierr,
			 (U_fp) F77_NAME(mmpy2), (U_fp) F77_NAME(smxpy2));
    } else if (*level == 4) {
	F77_CALL(blkfct)(m, nsuper, xsuper, snode, split,  xlindx, lindx,
			 xlnz, lnz, &iwsiz, iwork, &tmpsiz, tmpvec, ierr,
			 (U_fp) F77_NAME(mmpy4), (U_fp) F77_NAME(smxpy4));
    } else if (*level == 8) { /* the one used from R */
	F77_CALL(blkfct)(m, nsuper, xsuper, snode, split,  xlindx, lindx,
			 xlnz, lnz, &iwsiz, iwork, &tmpsiz, tmpvec, ierr,
			 (U_fp) F77_NAME(mmpy8), (U_fp) F77_NAME(smxpy8));
    }
    timewd[6] += F77_CALL(gtimer)() - timbeg;

    if (*ierr == -1) { *ierr = 10; return; }
    if (*ierr == -2) { *ierr = 11; return; }
    if (*ierr == -3) { *ierr = 12; return; }
    /* this is new [ from blkfct() -> blkfc2() in ./cholesky.c ] : */
    if (*ierr == -17) { *ierr = 17; return; }

} /* chlfct_ */

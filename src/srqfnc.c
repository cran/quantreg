/*-*- mode: C; kept-old-versions: 12;  kept-new-versions: 20; -*-
 *
 * srqfnc.f -- translated by f2c (version 20031025) and by
 * $Id: f2c-clean,v 1.10 2002/03/28 16:37:27 maechler Exp $
 * plus extended manual code cleaning by Martin Maechler, ETH Zurich
 */

#include <Rmath.h>

/* Table of constant values */

static int c__1 = 1;
static double c_9995 = .9995;


/* BLAS : */
#include <R_ext/BLAS.h>
/* double ddot_(int *, double *, int *, double *, int *); */
/* int   daxpy_(int *, double*, double *, int *, double *, int *); */

#define DDOT(_n_, _X_, _Y_)  F77_CALL(ddot)(_n_, _X_, &c__1, _Y_, &c__1)

/* SparseM -- ./sparskit2.f : */
#include "sparseM.h"

/* Cholesky related : */
#include "cholesky.h"

/* advance declaration : */
static void
slpfnc(int *n1, int *m,
       int *nnza1, double *a1, int *ja1, int *ia1,
       double *ao1, int *jao1, int *iao1,
       int *n2,
       int *nnza2, double *a2, int *ja2, int *ia2,
       double *ao2, int *jao2, int *iao2,
       int *nnzdmax, double *d, int *jd, int *id,
       double *dsub, int *jdsub, int *nsubmax, int *lindx,
       int *xlindx, int *nnzlmax, double *lnz, int *xlnz,
       int *invp, int *perm, int *iwmax, int *iwork, int *colcnt,
       int *snode, int *xsuper, int *split,
       int *tmpmax, double *tmpvec, double *rhs, double *newrhs,
       int *cachsz, int *level, double *x1, double *x2,
       double *s, double *u, double *c1, double *c2,
       double *y, double *b, double *r2, double *z1,
       double *z2, double *w, double *q1, double *q2,
       int *nnzemax, double *e, int *je, int *ie,
       int *nnzgmax, double *g, int *jg, int *ig,
       int *nnzhmax, double *h, int *jh, int *ih,
       double *dy, double *dx1, double *dx2, double *ds,
       double *dz1, double * dz2, double *dw,
       double *dxdz1, double *dxdz2, double *dsdw,
       double *xi1, double *xi2, int *maxn1n2,
       double *ww1, double *ww2, double *ww3,
       double *small, int *ierr, int *maxit, double *timewd);

static void
boundc(double *x1, double *dx1, double *x2,
       double *dx2, double *s, double *ds, double *z1,
       double *dz1, double *z2, double *dz2, double *w,
       double *dw, int *n1, int *n2, double *beta,
       double *deltap, double *deltad);

/* This is called from R : */
int
F77_SUB(srqfnc)(int *n1, int *m, int *nnza1,
	double *a1, int *ja1, int *ia1,
	double *ao1, int *jao1, int *iao1,
	int *n2, int *nnza2,
	double *a2, int *ja2, int *ia2,
	double *ao2, int *jao2, int *iao2,
	int *nnzdmax, double *d, int *jd, int *id,
	double *dsub, int *jdsub, int *nnzemax, double *e,
	int *je, int *ie, int *nnzgmax,
	double *g, int *jg, int *ig, int *nnzhmax,
	double *h, int *jh, int *ih,
	int *nsubmax, int *lindx, int *xlindx,
	int *nnzlmax, double *lnz, int *xlnz, int *iw,
	int *iwmax, int *iwork, int *xsuper, int *tmpmax,
	double *tmpvec, int *maxn1n2,
	double *ww1, double *wwm, double *wwn1, double *wwn2,
	int *cachsz, int *level, double *x1, double *x2,
	double *s, double *u,
	double *c1, double *c2, double *y, double *small,
	int *ierr, int *maxit, double *timewd)
{
    /* System generated locals */
    int iw_dim1, wwm_dim1, wwn1_dim1, wwn2_dim1;

    /* Parameter adjustments */
    wwm_dim1  = *m;	wwm -= wwm_dim1;
    wwn1_dim1 = *n1;	wwn1 -= wwn1_dim1;
    wwn2_dim1 = *n2;	wwn2 -= wwn2_dim1;
    iw_dim1   = *m;	iw -= iw_dim1;

    /* Function Body */
    slpfnc(n1, m, nnza1, a1, ja1, ia1, ao1, jao1, iao1,
	    n2, nnza2, a2, ja2, ia2, ao2, jao2, iao2, nnzdmax,
	    d, jd, id, dsub, jdsub,
	    nsubmax, lindx, xlindx, nnzlmax,
	    lnz, xlnz, &iw[iw_dim1], &iw[(iw_dim1 << 1)],
	    iwmax, iwork, &iw[iw_dim1 * 3], &iw[(iw_dim1 << 2)],
	    xsuper, &iw[iw_dim1 * 5], tmpmax, tmpvec,
	    &wwm[(wwm_dim1 << 1)], &wwm[wwm_dim1 * 3], cachsz, level,
	    x1, x2, s, u,
	    c1, c2, y, &wwm[wwm_dim1], &wwn2[wwn2_dim1],
	    &wwn1[wwn1_dim1], &wwn2[(wwn2_dim1 << 1)],
	    &wwn1[(wwn1_dim1 << 1)], &wwn1[wwn1_dim1 * 3],
	    &wwn2[wwn2_dim1 * 3], nnzemax, e, je, ie,
	    nnzgmax, g, jg, ig, nnzhmax, h, jh, ih,
	    &wwm[(wwm_dim1 << 2)],
	    &wwn1[(wwn1_dim1 << 2)], &wwn2[(wwn2_dim1 << 2)],
	    &wwn1[wwn1_dim1 * 5], &wwn1[wwn1_dim1 * 6],
	    &wwn2[wwn2_dim1 * 5],
	    &wwn1[wwn1_dim1 * 7], &wwn1[(wwn1_dim1 << 3)],
	    &wwn2[wwn2_dim1 * 6], &wwn1[wwn1_dim1 * 9],
	    &wwn1[wwn1_dim1 * 10], &wwn2[wwn2_dim1 * 7], maxn1n2, ww1,
	    &wwm[wwm_dim1 * 5], &wwm[wwm_dim1 * 6], small, ierr, maxit, timewd);
    return 0;
} /* srqfnc_ */

static void
slpfnc(int *n1, int *m,
       int *nnza1, double *a1, int *ja1, int *ia1,
       double *ao1, int *jao1, int *iao1,
       int *n2,
       int *nnza2, double *a2, int *ja2, int *ia2,
       double *ao2, int *jao2, int *iao2,
       int *nnzdmax, double *d, int *jd, int *id,
       double *dsub, int *jdsub, int *nsubmax, int *lindx,
       int *xlindx, int *nnzlmax, double *lnz, int *xlnz,
       int *invp, int *perm, int *iwmax, int *iwork, int *colcnt,
       int *snode, int *xsuper, int *split,
       int *tmpmax, double *tmpvec, double *rhs, double *newrhs,
       int *cachsz, int *level, double *x1, double *x2,
       double *s, double *u, double *c1, double *c2,
       double *y, double *b, double *r2, double *z1,
       double *z2, double *w, double *q1, double *q2,
       int *nnzemax, double *e, int *je, int *ie,
       int *nnzgmax, double *g, int *jg, int *ig,
       int *nnzhmax, double *h, int *jh, int *ih,
       double *dy, double *dx1, double *dx2, double *ds,
       double *dz1, double * dz2, double *dw,
       double *dxdz1, double *dxdz2, double *dsdw,
       double *xi1, double *xi2, int *maxn1n2,
       double *ww1, double *ww2, double *ww3,
       double *small, int *ierr, int *maxit, double *timewd)
{

/* Sparse implentation of LMS's interior point method via
    Ng-Peyton's sparse Cholesky factorization for sparse
    symmetric positive definite

 INPUT:
     n1 -- the number of row in the coefficient matrix A1'
     m -- the number of column in the coefficient matrix A1'
     nnza1 -- the number of non-zero elements in A'
     a1 -- an nnza1-vector of non-zero values of the design
          matrix (A1') stored in csr format
     ja1 -- an nnza1-vector of indices of the non-zero elements of
           the coefficient matrix
     ia1 -- an (n1+1)-vector of pointers to the begining of each
           row in a1 and ja1
     ao1 -- an nnza1-vector of work space for the transpose of
           the design matrix stored in csr format or the
           design matrix stored in csc format
     jao1 -- an nnza1-vector of work space for the indices of the
            transpose of the design matrix
     iao1 -- an (n1+1)-vector of pointers to the begining of each
            column in ao1 and jao1
     n2 -- the number of row in the constraint matrix A2'
     nnza2 -- the number of non-zero elements in A2'
     a2 -- an nnza2-vector of non-zero values of the contraint
          matrix (A2') stored in csr format
     ja2 -- an nnza2-vector of indices of the non-zero elements of
           the constraint matrix
     ia2 -- an (n2+1)-vector of pointers to the begining of each
           row in a2 and ja2
     ao2 -- an nnza2-vector of work space for the transpose of
           the constraint matrix stored in csr format or the
           constraint matrix stored in csc format
     jao2 -- an nnza2-vector of work space for the indices of the
            transpose of the constraint matrix
     iao2 -- an (n2+1)-vector of pointers to the begining of each
            column in ao2 and jao2
     nnzdmax -- upper bound of the non-zero elements in A1A1'
     d -- an nnzdmax-vector of non-zero values used to store
          the transpose of the design matrix multiplied by the design
          matrix (A1A1') stored in csr format;
          also used to store A1Q1^(-1) and A2Q2^(-1) later
     jd -- an nnzdmax-vector of indices in d
     id -- an (m+1)-vector of pointers to the begining of each
           row in d and jd
     dsub -- the values of d excluding the diagonal elements
     jdsub -- the indices to dsub
     nsubmax -- upper bound of the dimension of lindx
     lindx -- an nsub-vector of integer which contains, in
           column major order, the row subscripts of the nonzero
           entries in L in a compressed storage format
     xlindx -- an (m+1)-vector of int of pointers for lindx
     nnzlmax -- the upper bound of the non-zero entries in
                L stored in lnz, including the diagonal entries
     lnz -- First contains the non-zero entries of d; later
            contains the entries of the Cholesky factor
     xlnz -- column pointer for L stored in lnz
     invp -- an n1-vector of int of inverse permutation
             vector
     perm -- an n1-vector of int of permutation vector
     iw -- int work array of length m
     iwmax -- upper bound of the general purpose int
              working storage iwork; set at 7*m+3
     iwork -- an iwsiz-vector of int as work space
     colcnt -- array of length m, containing the number of
               non-zeros in each column of the factor, including
               the diagonal entries
     snode -- array of length m for recording supernode
              membership
     xsuper -- array of length m+1 containing the supernode
               partitioning
     split -- an m-vector with splitting of supernodes so that
              they fit into cache
     tmpmax -- upper bound of the dimension of tmpvec
     tmpvec -- a tmpmax-vector of temporary vector
     rhs -- m-vector to store the rhs
     newrhs -- extra work vector for right-hand side and
               solution
     cachsz -- size of the cache (in kilobytes) on the target
               machine
     level -- level of loop unrolling while performing numerical
              factorization
     x1 -- an n1-vector, the initial feasible solution for the primal
           solution that corresponds to the design matrix A1'
     x2 -- an n2-vector, the initial feasible solution for the primal
           solution that corresponds to the constraint matrix A2'
     s -- an n1-vector
     u -- an n1-vector of the upper bound for x1
     c1 -- an n1-vector in the primal; negative response in the
           regression quantile setting
     c2 -- an n2-vector, the negative rhs of the inequality constraint
     y -- an m-vector, the initial dual solution
     b -- an n1-vector, usualy the rhs of the equality constraint
          X'a = (1-tau)X'e in the rq setting
     r2 -- an n2-vector of residuals
     z1 -- an n1-vector of the dual slack variable
     z2 -- an n2-vector
     w -- an n-vector
     q1 -- an n1-vector of work array containing the diagonal
          elements of the Q1^(-1) matrix
     q2 -- an n2-vector of work array containing the diagonal
          elements of the Q2^(-1) matrix
     e -- an nnzdmax-vector containing the non-zero entries of
          A1Q1^(-1)A1' stored in csr format
     je -- an nnzdmax-vector of indices for e
     ie -- an (m+1)-vector of pointers to the begining of each
           row in e and je
     nnzgmax -- upper bound of the non-zero elements in g,jg
     g -- an nnzgmax-vector containing the non-zero entries of
          A2Q2^(-1)A2' stored in csr format
     jg -- an nnzgmax-vector of indices for g
     ig -- an (m+1)-vector of pointers to the begining of each
           row in g and jg
     nnzhmax -- upper bound of the non-zero elements in h,jh
     h -- an nnzhmax-vector containing the non-zero entries of
          AQ^(-1)A' stored in csr format
     jh -- an nnzhmax-vector of indices for h
     ih -- an (m+1)-vector of pointers to the begining of each
           row in h and jh
     dy -- an m-vector of work array
     dx1 -- an n1-vector of work array
     dx2 -- an n2-vector of work array
     ds -- an n1-vector of work array
     dz1 -- an n1-vector of work array
     dz2 -- an n2-vector of work array
     dw -- an n1-vector of work array
     dxdz1 -- an n1-vector of work array
     dxdz2 -- an n2-vector of work array
     dsdw -- an n1-vector of work arry
     xi1 -- an n1-vector of work array
     xi2 -- an n2-vector of work array
     xinv1 -- an n1-vector of work array
     xinv2 -- an n2-vector of work array
     sinv -- work array
     maxn1n2 -- max(n1,n2)
     ww1 -- an maxn1n2-vector of work array
     ww2 -- an m-vector of work array
     ww3 -- an m-vector of work array
     small -- convergence tolerance for interior algorithm
     ierr -- error flag :
       1 -- insufficient storage (work space) when calling extract;
       3 -- insufficient storage in iwork when calling ordmmd;
       4 -- insufficient storage in iwork when calling sfinit;
       5 -- nnzl > nnzlmax when calling sfinit
       6 -- nsub > nsubmax when calling sfinit
       7 -- insufficient work space in iwork when calling symfct
       8 -- inconsistancy in input when calling symfct
       9 -- tmpsiz > tmpmax when calling symfct; increase tmpmax
       10 -- nonpositive diagonal encountered    when calling blkfct
       11 -- insufficient work storage in tmpvec when calling blkfct
       12 -- insufficient work storage in iwork  when calling blkfct
       13 -- nnzd > nnzdmax in e,je when calling amub
       14 -- nnzd > nnzdmax in g,jg when calling amub
       15 -- nnzd > nnzdmax in h,jh when calling aplb
       17 -- tiny diagonals replaced with Inf when calling blkfct()
       		{new error code, was confounded with '10'}

     maxit -- maximal number of iterations; on return: the number of iterations

     timewd -- vector of length 7: [7]: amount of time for this subroutine
               [1:6] time info for the phases of cholfct() only.

 OUTPUT:
     y -- an m-vector of primal solution
*/

    /* Local variables */
    int i, it, nnzd, nsuper, nnzdsub;
    double mu, gap, deltad, deltap;
    double timbeg, timend;

/* Parameter adjustments */
    --dsdw;
    --dxdz1;
    --dxdz2;
    --dz1; --dz2;
    --dx1; --dx2;
    --ds;
    --dw;
    --dy;
    --c1;  --c2;
    --q1;  --q2;
    --x1;  --x2;
    --z1;  --z2;
    --r2;
    --xi1; --xi2;
    --b;
    --s;
    --u;
    --w;
    --y;

    --ww1;
    --ww2;
    --ww3;

    --newrhs;
    --rhs;

    --perm;
    --invp;

    --dsub;
    --timewd;

    /* Function Body */
    for (i = 1; i <= 7; ++i)
	timewd[i] = 0.;

/*  Compute the initial gap */

    gap =
	DDOT(n1, &z1[1], &x1[1]) +
	DDOT(n2, &z2[1], &x2[1]) +
	DDOT(n1, &w[1],  &s[1]);

/*  Start iteration */
    it = 0;
    while(gap >= *small && it <= *maxit) {

	++it;

/*  Create the diagonal matrix Q1^(-1) stored in q1 as an n1-vector,
    the diagonal matrix Q2^(-1) stored in q2 as an n2-vector,
    and store the residuals in r1 in ds, and r3 in dy temporarily,
    and r2 in r2 permanently
*/

	/* amux: obtain A1x1 and store the value in ww2 */
	F77_CALL(amux)(m, &x1[1], &ww2[1], ao1, jao1, iao1);

	/* amux: obtain A2x2 and store the value in ww3 */
	F77_CALL(amux)(m, &x2[1], &ww3[1], ao2, jao2, iao2);

	/*    obtain  A2'y and store it temporarily in r2 */
	F77_CALL(amux)(n2, &y[1], &r2[1], a2, ja2, ia2);

	for (i = 1; i <= *n1; ++i) {
	    q1[i] = 1. / (z1[i] / x1[i] + w[i] / s[i]);
	    ds[i] = z1[i] - w[i];
	}
	for (i = 1; i <= *n2; ++i) {
	    q2[i] = x2[i] / z2[i];
	    r2[i] = c2[i] - r2[i];
	}
	for (i = 1; i <= *m; ++i) {
	    dy[i] = b[i] - ww2[i] - ww3[i];
	}

/*  Obtain AQA = A1Q1^(-1)A1' + A2Q2^(-1)A2' in 5 steps :

 *  Step1: Obtain A1Q1^(-1) and store the values in d,jd,id in csr format
 *         Also compute A1Q1^(-1)r1 and store the values in ww2 to be used
 *         to generate r3;
 *  Step2: Compute A1Q1^(-1)A1' and store the values in e,je,ie
 *  Step3: Obtain A2Q2^(-1) and store the values in d,jd,id in csr format
 *         Also compute A2Q2^(-1)r2 and store the values in in ww3 to
 *         be used to generate r3;
 *  Step4: Compute A2Q2^(-1)A2' and store the value in g,jg,ig
 *  Step5: Compute AQA and store the values in h,jh,ih
 */

/*   Step 1 */

	F77_CALL(amudia)(m, &c__1, ao1, jao1, iao1, &q1[1], d, jd, id);
	F77_CALL(amux)(m, &ds[1], &ww2[1], d, jd, id);

/*   Step 2 */

	F77_CALL(amub)(m, m, &c__1, d, jd, id, a1, ja1, ia1,
		       e, je, ie, nnzemax, iwork, ierr);
	if (*ierr) { *ierr = 13; goto L100; }

/*   Step 3 */

	F77_CALL(amudia)(m, &c__1, ao2, jao2, iao2, &q2[1], d, jd, id);
	F77_CALL(amux)(m, &r2[1], &ww3[1], d, jd, id);

/*   Step 4 */

	F77_CALL(amub)(m, m, &c__1, d, jd, id, a2, ja2, ia2,
		       g, jg, ig, nnzgmax, iwork, ierr);
	if (*ierr) { *ierr = 14; goto L100; }

/*   Step 5 */

	F77_CALL(aplb)(m, m, &c__1, e, je, ie, g, jg, ig,
		       h, jh, ih, nnzhmax, iwork, ierr);
	if (*ierr) { *ierr = 15; goto L100; }

	/* Generate rhs = r3 + A1Q1^(-1) r1 + A2Q2^(-1) r2 : */

	for (i = 1; i <= *m; ++i) {
	    rhs[i] = dy[i] + ww2[i] + ww3[i];
	}

/*  Extract the non-diagonal structure of h,jh,ih and store in dsub,jdsub */

	nnzd = ih[*m] - 1;
	nnzdsub = nnzd - *m;
	i = *nnzhmax + 1;
	F77_CALL(extract)(h, jh, ih, &dsub[1], jdsub, m, nnzhmax, &i, ierr);

	if (*ierr) { *ierr = 1; goto L100; }

/*  Compute dy = (AQ^(-1)A')^(-1)rhs; result returned via dy

Call chlfct to perform Cholesky's decomposition of h,jh,ih */

	F77_CALL(chlfct)(m, xlindx, lindx, &invp[1], &perm[1], iwork, &nnzdsub,
			 jdsub, colcnt, &nsuper, snode, xsuper, nnzlmax,
			 nsubmax, xlnz, lnz, ih, jh, h, cachsz,
			 tmpmax, level, tmpvec, split, ierr, &it, &timewd[1]);

	if (*ierr) { goto L100; }

/* Call blkslv: Numerical solution for the new rhs stored in rhs */

	for (i = 1; i <= *m; ++i) {
	    newrhs[i] = rhs[perm[i]];
	}
	timbeg = F77_CALL(gtimer)();
	F77_CALL(blkslv)(&nsuper, xsuper, xlindx, lindx, xlnz, lnz, &newrhs[1]);
	timend = F77_CALL(gtimer)();
	timewd[7] = timewd[7] + timend - timbeg;
	for (i = 1; i <= *m; ++i) {
	    dy[i] = newrhs[invp[i]];
	}

/*  Compute dx1 = Q1^(-1)(A1'dy - r1), ds = -dx1, dz1, dz2  and dw */

	F77_CALL(amux)(n1, &dy[1], &dx1[1], a1, ja1, ia1);
	F77_CALL(amux)(n2, &dy[1], &dx2[1], a2, ja2, ia2);
	for (i = 1; i <= *n1; ++i) {
	    dx1[i] = q1[i] * (dx1[i] - ds[i]);
	    ds[i] = -dx1[i];
	    dz1[i] = -z1[i] * (dx1[i] / x1[i] + 1.);
	    dw[i] = -w[i] * (ds[i] / s[i] + 1.);
	}
	for (i = 1; i <= *n2; ++i) {
	    dx2[i] = q2[i] * (dx2[i] - r2[i]);
	    dz2[i] = -z2[i] * (dx2[i] / x2[i] + 1.);
	}

/* Compute the maximum allowable step lengths */

	boundc(&x1[1], &dx1[1], &x2[1], &dx2[1],
	       &s[1],  &ds[1],  &z1[1], &dz1[1],
	       &z2[1], &dz2[1], &w[1], &dw[1],
	       n1, n2, &c_9995, &deltap, &deltad);

	if (deltap * deltad < 1.) {

	    /* Update mu */
	    double g1;

	    mu = DDOT(n1, &z1[1], &x1[1]) +
		DDOT(n2, &z2[1], &x2[1]) +
		DDOT(n1, &w[1], &s[1]);

	    g1 = mu +
		deltap *          DDOT(n1, &z1[1], &dx1[1]) +
		deltad *          DDOT(n1, &dz1[1], &x1[1]) +
		deltad * deltap * DDOT(n1, &dz1[1], &dx1[1]) +
		deltap *          DDOT(n2, &z2[1], &dx2[1]) +
		deltad *          DDOT(n2, &dz2[1], &x2[1]) +
		deltad * deltap * DDOT(n2, &dz2[1], &dx2[1]) +
		deltap *          DDOT(n1, &w[1], &ds[1]) +
		deltad *          DDOT(n1, &dw[1], &s[1]) +
		deltad * deltap * DDOT(n1, &dw[1], &ds[1]);

	    g1 /= mu;
	    mu *=  g1 * g1 * g1 / (2. * *n1 + *n2);

/* Compute dx1dz1, dx2dz2 and dsdw */

	    for (i = 1; i <= *n1; ++i) {
		dxdz1[i] = dx1[i] * dz1[i];
		dsdw[i] = ds[i] * dw[i];
		xi1[i] = dxdz1[i] / x1[i] - dsdw[i] / s[i] -
		     mu * (    1. / x1[i] -      1. / s[i]);
		ww1[i] = q1[i] * xi1[i];
	    }

/* Compute A1Q1^(-1)(X1^(-1)*dx1dz1 - S^(-1)*dsdw - mu(X1^(-1) - S^(-1))) and
   store it in ww2 temporarily */

	    F77_CALL(amux)(m, &ww1[1], &ww2[1], ao1, jao1, iao1);
	    for (i = 1; i <= *n2; ++i) {
		dxdz2[i] = dx2[i] * dz2[i];
		xi2[i] = (dxdz2[i] - mu) / x2[i];
		ww1[i] = q2[i] * xi2[i];
	    }

/* Compute A2Q2^(-1)(X2^(-1)*dx2dz2 - mu X2^(-1)) and store it in ww3 temporarily */

	    F77_CALL(amux)(m, &ww1[1], &ww3[1], ao2, jao2, iao2);
	    for (i = 1; i <= *m; ++i) {
		rhs[i] +=  ww2[i] + ww3[i];
	    }


/* Compute (AQ^(-1)A')^(-1)rhs and return the result in dy

Call blkslv: Numerical solution for the new rhs stored in rhs */

	    for (i = 1; i <= *m; ++i) {
		newrhs[i] = rhs[perm[i]];
	    }
	    timbeg = F77_CALL(gtimer)();
	    F77_CALL(blkslv)(&nsuper, xsuper, xlindx, lindx, xlnz, lnz, &newrhs[1]);
	    timend = F77_CALL(gtimer)();
	    timewd[7] = timewd[7] + timend - timbeg;
	    for (i = 1; i <= *m; ++i) {
		dy[i] = newrhs[invp[i]];
	    }

/*  Compute dx1=Q1^(-1)(A1'dy-X1^(-1)*dx1dz1-S^(-1)*dsdw
    -mu*(X1^(-1)-S^(-1))-r1), ds = -dx1, dz1, dz2  and dw */

	    F77_CALL(amux)(n1, &dy[1], &dx1[1], a1, ja1, ia1);
	    F77_CALL(amux)(n2, &dy[1], &dx2[1], a2, ja2, ia2);
	    for (i = 1; i <= *n1; ++i) {
		dx1[i] = q1[i] * (dx1[i] - xi1[i] - z1[i] + w[i]);
		ds[i] = -dx1[i];
		dz1[i] = -z1[i] + (mu - z1[i] * dx1[i] - dxdz1[i]) / x1[i];
		dw[i] = -w[i] + (mu - w[i] * ds[i] - dsdw[i]) / s[i];
	    }
	    for (i = 1; i <= *n2; ++i) {
		dx2[i] = q2[i] * (dx2[i] - xi2[i] - r2[i]);
		dz2[i] = -z2[i] + (mu - z2[i] * dx2[i] - dxdz2[i]) / x2[i];
	    }

/* Compute the maximum allowable step lengths */

	    boundc(&x1[1], &dx1[1], &x2[1], &dx2[1], &s[1], &ds[1],
		   &z1[1], &dz1[1], &z2[1], &dz2[1], &w[1], &dw[1],
		   n1, n2, &c_9995, &deltap, &deltad);
	}

/* Take the step */

	F77_CALL(daxpy)(n1, &deltap, &dx1[1], &c__1, &x1[1], &c__1);
	F77_CALL(daxpy)(n2, &deltap, &dx2[1], &c__1, &x2[1], &c__1);
	F77_CALL(daxpy)(n1, &deltap, &ds[1],  &c__1, &s[1],  &c__1);
	F77_CALL(daxpy)(n1, &deltad, &dw[1],  &c__1, &w[1],  &c__1);
	F77_CALL(daxpy)(n1, &deltad, &dz1[1], &c__1, &z1[1], &c__1);
	F77_CALL(daxpy)(n2, &deltad, &dz2[1], &c__1, &z2[1], &c__1);
	F77_CALL(daxpy)(m,  &deltad, &dy[1],  &c__1, &y[1],  &c__1);
	gap =
	    DDOT(n1, &z1[1], &x1[1]) +
	    DDOT(n2, &z2[1], &x2[1]) +
	    DDOT(n1, &w[1],  &s[1] );

    } /* end {while} iterations */

L100:
    *maxit = it;
    return;
} /* slpfnc */

static void
boundc(double *x1, double *dx1, double *x2,
       double *dx2, double *s, double *ds, double *z1,
       double *dz1, double *z2, double *dz2, double *w,
       double *dw, int *n1, int *n2, double *beta,
       double *deltap, double *deltad)
{
    int i;

    *deltap = R_PosInf;
    *deltad = R_PosInf;

    for (i = 0; i < *n1; ++i) {
	if(dx1[i] < 0) *deltap = fmin2(*deltap, -x1[i] / dx1[i]);
	if(ds[i]  < 0) *deltap = fmin2(*deltap, -s[i] / ds[i]);
	if(dz1[i] < 0) *deltad = fmin2(*deltad, -z1[i] / dz1[i]);
	if(dw[i]  < 0) *deltad = fmin2(*deltad, -w[i] / dw[i]);
    }
    for (i = 0; i < *n2; ++i) {
	if(dx2[i] < 0) *deltap = fmin2(*deltap, -x2[i] / dx2[i]);
	if(dz2[i] < 0) *deltad = fmin2(*deltad, -z2[i] / dz2[i]);
    }
    *deltap = fmin2(1., *beta * *deltap);
    *deltad = fmin2(1., *beta * *deltad);
    return;
} /* boundc */

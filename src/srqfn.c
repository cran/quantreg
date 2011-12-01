/*-*- mode: C; kept-old-versions: 12;  kept-new-versions: 20; -*-
 *
 * srqfn.f -- translated by f2c (version 20031025) and by
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
slpfn(int *n, int *m, int *nnza, double *a, int *ja, int *ia,
      double *ao, int *jao, int *iao, int *nnzdmax, double *d, int *jd, int *id,
      double *dsub, int *jdsub, int *nsubmax, int *lindx, int *xlindx,
      int *nnzlmax, double *lnz, int *xlnz,
      int *invp, int *perm, int *iwmax, int *iwork,
      int *colcnt, int *snode, int *xsuper, int *split,
      int *tmpmax, double *tmpvec, double *newrhs, int *cachsz,
      int *level, double *x, double *s, double *u,
      double *c, double *y, double *b, double *r,
      double *z, double *w, double *q, int *nnzemax,
      double *e, int *je, int *ie,
      double *dy, double *dx, double *ds, double *dz, double *dw,
      double *dxdz, double *dsdw, double *xi, double *xinv, double *sinv,
      double *ww1, double *ww2, double *small, int *ierr,
      int *maxit, double *timewd);

static void
bound(double *x, double *dx, double *s,
      double *ds, double *z, double *dz, double *w,
      double *dw, int *n, double *beta, double *deltap,
      double *deltad);


/*  called from R's rq.fit.sfn()  in ../R/sfn.R */
int
srqfn_(int *n, int *m, int *nnza,
       double *a, int *ja, int *ia,
       double *ao, int *jao, int *iao,
       int *nnzdmax, double *d, int *jd, int *id,
	double *dsub, int *jdsub,
       int *nnzemax, double *e, int *je, int *ie,
       int *nsubmax, int *lindx, int *xlindx,
       int *nnzlmax, double *lnz, int *xlnz, int *iw,
       int *iwmax, int *iwork, int *xsuper, int *tmpmax,
       double *tmpvec, double *wwm, double *wwn, int *cachsz,
       int *level, double *x, double *s, double *u,
       double *c, double *y, double *b, double *small,
       int *ierr, int *maxit, double *timewd)
{
    /* System generated locals */
    int iw_dim1, wwm_dim1, wwn_dim1;

    /* Parameter adjustments */
    wwn_dim1 = *n;    wwn -= wwn_dim1;
    wwm_dim1 = *m;    wwm -= wwm_dim1;
    iw_dim1 = *m;     iw -= iw_dim1;


    /* Function Body */
    slpfn(n, m, nnza, a, ja, ia, ao, jao, iao,
	  nnzdmax, d, jd, id, dsub, jdsub, nsubmax,
	  lindx, xlindx, nnzlmax, lnz, xlnz,
	  &iw[iw_dim1], &iw[(iw_dim1 << 1)],
	  iwmax, iwork,
	  &iw[iw_dim1 * 3], &iw[(iw_dim1 << 2)], xsuper, &iw[iw_dim1 * 5],
	  tmpmax, tmpvec, &wwm[(wwm_dim1 << 1)], cachsz, level,
	  x, s, u, c, y, b,
	  &wwn[wwn_dim1    ], &wwn[(wwn_dim1 << 1)],
	  &wwn[wwn_dim1 * 3], &wwn[(wwn_dim1 << 2)],
	  nnzemax, e, je, ie,
	  &wwm[wwm_dim1 * 3], &wwn[wwn_dim1 * 5],
	  &wwn[wwn_dim1 * 6], &wwn[wwn_dim1 * 7],
	  &wwn[(wwn_dim1 << 3)], &wwn[wwn_dim1 * 9],
	  &wwn[wwn_dim1 * 10], &wwn[wwn_dim1 * 11],
	  &wwn[wwn_dim1 * 12], &wwn[wwn_dim1 * 13],
	  &wwn[wwn_dim1 * 14],
	  &wwm[wwm_dim1], small, ierr, maxit, timewd);
    return 0;
} /* srqfn_ */

static void
slpfn(int *n, int *m, int *nnza, double *a, int *ja, int *ia,
      double *ao, int *jao, int *iao, int *nnzdmax, double *d, int *jd, int *id,
      double *dsub, int *jdsub, int *nsubmax, int *lindx, int *xlindx,
      int *nnzlmax, double *lnz, int *xlnz,
      int *invp, int *perm, int *iwmax, int *iwork,
      int *colcnt, int *snode, int *xsuper, int *split,
      int *tmpmax, double *tmpvec, double *newrhs, int *cachsz,
      int *level, double *x, double *s, double *u,
      double *c, double *y, double *b, double *r,
      double *z, double *w, double *q, int *nnzemax,
      double *e, int *je, int *ie,
      double *dy, double *dx, double *ds, double *dz, double *dw,
      double *dxdz, double *dsdw, double *xi, double *xinv, double *sinv,
      double *ww1, double *ww2, double *small, int *ierr,
      int *maxit, double *timewd)
{
/* Sparse implentation of LMS's interior point method via
    Ng-Peyton's sparse Cholesky factorization for sparse
    symmetric positive definite

 INPUT:
     n -- the number of rows in the coefficient matrix A'
     m -- the number of columns in the coefficient matrix A'
     nnza -- the number of non-zero elements in A'
     a -- an nnza-vector of non-zero values of the design
          matrix (A') stored in csr format
     ja -- an nnza-vector of indices of the non-zero elements of
           the coefficient matrix
     ia -- an (n+1)-vector of pointers to the begining of each
           row in a and ja
     ao -- an nnza-vector of work space for the transpose of
           the design matrix stored in csr format or the
           design matrix stored in csc format
     jao -- an nnza-vector of work space for the indices of the
            transpose of the design matrix
     iao -- an (n+1)-vector of pointers to the begining of each
            column in ao and jao
     nnzdmax -- upper bound of the non-zero elements in AA'
     d -- an nnzdmax-vector of non-zero values of AQ^(-1)
     jd -- an nnzdmax-vector of indices in d
     id -- an (m+1)-vector of pointers to the begining of each
           row in d and jd
     dsub -- the values of e excluding the diagonal elements
     jdsub -- the indices to dsub
     nsubmax -- upper bound of the dimension of lindx
     lindx -- an nsub-vector of integer which contains, in
           column major oder, the row subscripts of the nonzero
           entries in L in a compressed storage format
     xlindx -- an (m+1)-vector of int of pointers for lindx
     nnzlmax -- the upper bound of the non-zero entries in
                L stored in lnz, including the diagonal entries
     lnz -- First contains the non-zero entries of d; later
            contains the entries of the Cholesky factor
     xlnz -- column pointer for L stored in lnz
     invp -- an n-vector of int of inverse permutation vector
     perm -- an n-vector of int of permutation vector
     iw -- int work array of length m
     iwmax -- upper bound of the general purpose int
              working storage iwork; set at 7*m+3
     iwork -- an iwsiz-vector of int as work space
     colcnt -- array of length m, containing the number of
               non-zeros in each column of the factor, including
               the diagonal entries
     snode -- array of length m for recording supernode membership
     xsuper -- array of length m+1 containing the supernode
               partitioning
     split -- an m-vector with splitting of supernodes so that
              they fit into cache
     tmpmax -- upper bound of the dimension of tmpvec
     tmpvec -- a tmpmax-vector of temporary vector
     newrhs -- extra work vector for right-hand side and solution
     cachsz -- size of the cache (in kilobytes) on the target machine
     level -- level of loop unrolling while performing numerical
              factorization
     x -- an n-vector, the initial feasible solution in the primal
              that corresponds to the design matrix A'
     s -- an n-vector
     u -- an n-vector of upper bound for x
     c -- an n-vector, usually the "negative" of
          the pseudo response
     y -- an m-vector, the initial dual solution
     b -- an n-vector, usualy the rhs of the equality constraint
          X'a = (1-tau)X'e in the rq setting
     r -- an n-vector of residuals
     z -- an n-vector of the dual slack variable
     w -- an n-vector
     q -- an n-vector of work array containing the diagonal
          elements of the Q^(-1) matrix
     nnzemax -- upper bound of the non-zero elements in AA'
     e -- an nnzdmax-vector containing the non-zero entries of
          AQ^(-1)A' stored in csr format
     je -- an nnzemax-vector of indices for e
     ie -- an (m+1)-vector of pointers to the begining of each row in e and je
     dy -- work array
     dx -- work array
     ds -- work array
     dz -- work array
     dw -- work array
     dxdz -- work array
     dsdw -- work arry
     xi -- work array
     xinv -- work array
     sinv -- work array
     ww1 -- work array
     ww2 -- work array
     small -- convergence tolerance for interior algorithm
     ierr -- error flag :
       1 -- insufficient storage (work space) when calling extract;
       2 -- nnzd > nnzdmax
       3 -- insufficient storage in iwork when calling ordmmd;
       4 -- insufficient storage in iwork when calling sfinit;
       5 -- nnzl > nnzlmax when calling sfinit
       6 -- nsub > nsubmax when calling sfinit
       7 -- insufficient work space in iwork when calling symfct
       8 -- inconsistancy in input when calling symfct
       9 -- tmpsiz > tmpmax when calling bfinit; increase tmpmax
       10 -- nonpositive diagonal encountered when calling blkfct(),
             the matrix is not positive definite
       11 -- insufficient work storage in tmpvec when calling blkfct()
       12 -- insufficient work storage in iwork  when calling blkfct()
       17 -- tiny diagonals replaced with Inf when calling blkfct()
       		{new error code, was confounded with '10'}

     maxit -- maximal number of iterations; on return: the number of iterations

     timewd -- vector of length 7: [7]: amount of time for this subroutine
               [1:6] time info for the phases of cholfct() only.

 OUTPUT:
     y -- an m-vector of primal solution
*/

    /* Local variables */
    int i, it, nnzd, nnzdsub, nsuper;
    double deltad, deltap, mu, gap;
    double timbeg;


/* Parameter adjustments */
    --sinv;
    --xinv;
    --xi;
    --dsdw;
    --dxdz;
    --dw;
    --dz;
    --ds;
    --dx;
    --q;
    --w;
    --z;
    --r;
    --c;
    --u;
    --s;
    --x;

    --ww1;
    --ww2;
    --dy;
    --b;
    --y;
    --newrhs;

    --perm;
    --invp;

    --dsub;
    --timewd;

    /* Function Body */
    for (i = 1; i <= 7; ++i)
	timewd[i] = 0.;

    nnzd = ie[*m] - 1;
    nnzdsub = nnzd - *m;

/*  Compute the initial gap */

    gap = DDOT(n, &z[1], &x[1]) + DDOT(n, &w[1], &s[1]);

/*  Start iteration */

    it = 0;
    while(gap >= *small && it <= *maxit) {

	++it;

/*  Create the diagonal matrix Q^(-1) stored in q as an n-vector
    and update the residuals in r */

	for (i = 1; i <= *n; ++i) {
	    q[i] = 1. / (z[i] / x[i] + w[i] / s[i]);
	    r[i] = z[i] - w[i];
	}

/*  Obtain AQ^(-1) and store in d,jd,id in csr format */

	F77_CALL(amudia)(m, &c__1, ao, jao, iao, &q[1], d, jd, id);

/*  Obtain AQ^(-1)A' and store in e,je,ie in csr format */

	F77_CALL(amub)(m, m, &c__1, d, jd, id,
		       a, ja, ia, e, je, ie, nnzemax, iwork, ierr);

	if (*ierr != 0) { *ierr = 2; goto L100; }

/*  Extract the non-diagonal structure of e,je,ie and store in dsub,jdsub */

	i = *nnzemax + 1;
	F77_CALL(extract)(e, je, ie, &dsub[1], jdsub, m, nnzemax, &i, ierr);

	if (*ierr != 0) { *ierr = 1; goto L100; }

/* Compute b - Ax + AQ^(-1)r and store it in c in two steps
   First: store Ax in ww2 */
	F77_CALL(amux)(m, &x[1], &ww2[1], ao, jao, iao);

/*  Second: save AQ^(-1)r in c temporarily */

	F77_CALL(amux)(m, &r[1], &c[1], d, jd, id);
	for (i = 1; i <= *m; ++i) {
	    c[i] = b[i] - ww2[i] + c[i];
	}

/*  Compute dy = (AQ^(-1)A')^(-1)(b-Ax+AQ^(-1)r); result returned via dy */

	/* chlfct: perform Cholesky's decomposition of e,je,ie */

	F77_CALL(chlfct)(m, xlindx, lindx, &invp[1], &perm[1], iwork, &nnzdsub,
			 jdsub, colcnt, &nsuper, snode, xsuper, nnzlmax,
			 nsubmax, xlnz, lnz, ie, je, e, cachsz, tmpmax,
			 level, tmpvec, split, ierr, &it, &timewd[1]);

	if (*ierr != 0) { /* return the chlfct() error code in {3:12}: */ goto L100; }

	for (i = 1; i <= *m; ++i)
	    newrhs[i] = c[perm[i]];

	/* blkslv: Numerical solution for the new rhs stored in b */
	timbeg = F77_CALL(gtimer)();
	F77_CALL(blkslv)(&nsuper, xsuper, xlindx, lindx, xlnz, lnz, &newrhs[1]);
	timewd[7] += F77_CALL(gtimer)() - timbeg;
	for (i = 1; i <= *m; ++i)
	    dy[i] = newrhs[invp[i]];

/*  Compute dx = Q^(-1)(A'dy - r), ds = -dx, dz  and dw */

	F77_CALL(amux)(n, &dy[1], &dx[1], a, ja, ia);
	for (i = 1; i <= *n; ++i) {
	    dx[i] = q[i] * (dx[i] - r[i]);
	    ds[i] = -dx[i];
	    dz[i] = -z[i] * (dx[i] / x[i] + 1.);
	    dw[i] = -w[i] * (ds[i] / s[i] + 1.);
	}

/* Compute the maximum allowable step lengths */

	bound(&x[1], &dx[1], &s[1], &ds[1], &z[1], &dz[1], &w[1], &dw[1],
	      n, &c_9995, &deltap, &deltad);

	if (deltap * deltad < 1.) {

	    /* Update mu */
	    double
		g =		      DDOT(n, &z[1], &x[1]) +
		    deltap *	      DDOT(n, &z[1], &dx[1]) +
		    deltad *	      DDOT(n, &dz[1], &x[1]) +
		    deltad * deltap * DDOT(n, &dz[1], &dx[1]) +
				      DDOT(n, &w[1], &s[1]) +
		    deltap *	      DDOT(n, &w[1], &ds[1]) +
		    deltad *	      DDOT(n, &dw[1], &s[1]) +
		    deltad * deltap * DDOT(n, &dw[1], &ds[1]);

	    mu = DDOT(n, &z[1], &x[1]) + DDOT(n, &w[1], &s[1]);
	    g /= mu;
	    mu *= g * g * g / (2. * (*n));

	    /* Compute dxdz and dsdw */

	    for (i = 1; i <= *n; ++i) {
		dxdz[i] = dx[i] * dz[i];
		dsdw[i] = ds[i] * dw[i];
		xinv[i] = 1. / x[i];
		sinv[i] = 1. / s[i];
		xi[i] = xinv[i] * dxdz[i] - sinv[i] * dsdw[i] -
		    mu * (xinv[i] - sinv[i]);
		ww1[i] = q[i] * xi[i];
	    }

	    /* Compute AQ^(-1)(dxdz - dsdw - mu(X^(-1) - S^(-1))) and
	       store it in ww2 temporarily */

	    F77_CALL(amux)(m, &ww1[1], &ww2[1], ao, jao, iao);
	    for (i = 1; i <= *m; ++i) {
		c[i] += ww2[i];
	    }

	    /* Compute dy and return the result in dy */

	    for (i = 1; i <= *m; ++i) {
		newrhs[i] = c[perm[i]];
	    }

	    /* Call blkslv: Numerical solution for the new rhs stored in b */

	    timbeg = F77_CALL(gtimer)();
	    F77_CALL(blkslv)(&nsuper, xsuper, xlindx, lindx, xlnz, lnz, &newrhs[1]);
	    timewd[7] += F77_CALL(gtimer)() - timbeg;
	    for (i = 1; i <= *m; ++i) {
		dy[i] = newrhs[invp[i]];
	    }

	    /* Compute dx = Q^(-1)(A'dy - r + mu(X^(-1) - S^(-1)) -dxdz + dsdw),
	       ds = -dx, dz  and dw */

	    F77_CALL(amux)(n, &dy[1], &dx[1], a, ja, ia);
	    for (i = 1; i <= *n; ++i) {
		dx[i] = q[i] * (dx[i] - xi[i] - r[i]);
		ds[i] = -dx[i];
		dz[i] = -z[i] + xinv[i] * (mu - z[i] * dx[i] - dxdz[i]);
		dw[i] = -w[i] + sinv[i] * (mu - w[i] * ds[i] - dsdw[i]);
	    }

	    /* Compute the maximum allowable step lengths */

	    bound(&x[1], &dx[1], &s[1], &ds[1], &z[1], &dz[1], &w[1], &dw[1],
		  n, &c_9995, &deltap, &deltad);
	}

	/* Take the step */

	F77_CALL(daxpy)(n, &deltap, &dx[1], &c__1, &x[1], &c__1);
	F77_CALL(daxpy)(n, &deltap, &ds[1], &c__1, &s[1], &c__1);
	F77_CALL(daxpy)(n, &deltad, &dw[1], &c__1, &w[1], &c__1);
	F77_CALL(daxpy)(n, &deltad, &dz[1], &c__1, &z[1], &c__1);
	F77_CALL(daxpy)(m, &deltad, &dy[1], &c__1, &y[1], &c__1);

	gap = DDOT(n, &z[1], &x[1]) + DDOT(n, &w[1], &s[1]);

    } /* end {while} iterations */


L100:
    *maxit = it;
    return;
} /* slpfn */

static void
bound(double *x, double *dx, double *s,
      double *ds, double *z, double *dz, double *w,
      double *dw, int *n, double *beta, double *deltap,
      double *deltad)
{
    int i;

    *deltap = R_PosInf;
    *deltad = R_PosInf;
    for (i = 0; i < *n; ++i) {
	if (dx[i] < 0.) *deltap = fmin2(*deltap, -x[i] / dx[i]);
	if (ds[i] < 0.) *deltap = fmin2(*deltap, -s[i] / ds[i]);
	if (dz[i] < 0.) *deltad = fmin2(*deltad, -z[i] / dz[i]);
	if (dw[i] < 0.) *deltad = fmin2(*deltad, -w[i] / dw[i]);
    }

    *deltap = fmin2(1., *beta * *deltap);
    *deltad = fmin2(1., *beta * *deltad);
    return;
} /* bound */

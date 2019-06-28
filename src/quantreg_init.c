#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void bootnp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(brutpow)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(combin)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(crqf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(crqfnb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kuantiles)(void *, void *, void *, void *);
extern void F77_NAME(penalty)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(powell)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pwy)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rls)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rqbr)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rqfnb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rqfnc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rqs)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sakj)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(srqfn)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(srqfnc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(wxy)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xys)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bootnp", (DL_FUNC) &bootnp, 13},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"brutpow",   (DL_FUNC) &F77_NAME(brutpow),   14},
    {"combin",    (DL_FUNC) &F77_NAME(combin),     7},
    {"crqf",      (DL_FUNC) &F77_NAME(crqf),      26},
    {"crqfnb",    (DL_FUNC) &F77_NAME(crqfnb),    18},
    {"kuantiles", (DL_FUNC) &F77_NAME(kuantiles),  4},
    {"penalty",   (DL_FUNC) &F77_NAME(penalty),   14},
    {"powell",    (DL_FUNC) &F77_NAME(powell),    17},
    {"pwy",       (DL_FUNC) &F77_NAME(pwy),       16},
    {"rls",       (DL_FUNC) &F77_NAME(rls),        7},
    {"rqbr",      (DL_FUNC) &F77_NAME(rqbr),      27},
    {"rqfnb",     (DL_FUNC) &F77_NAME(rqfnb),     13},
    {"rqfnc",     (DL_FUNC) &F77_NAME(rqfnc),     18},
    {"rqs",       (DL_FUNC) &F77_NAME(rqs),       15},
    {"sakj",      (DL_FUNC) &F77_NAME(sakj),      13},
    {"srqfn",     (DL_FUNC) &F77_NAME(srqfn),     45},
    {"srqfnc",    (DL_FUNC) &F77_NAME(srqfnc),    65},
    {"wxy",       (DL_FUNC) &F77_NAME(wxy),       18},
    {"xys",       (DL_FUNC) &F77_NAME(xys),       19},
    {NULL, NULL, 0}
};

void R_init_quantreg(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

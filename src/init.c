#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(code1)(void *, void *, void *, void *, void*, void*);
extern void F77_NAME(code2)(void *, void *, void *, void*, void*);
extern void F77_NAME(code3)(void *, void *, void *, void*);

static const R_FortranMethodDef FortranEntries[] = {
    {"code1", (DL_FUNC) &F77_NAME(code1), 6},
    {"code2", (DL_FUNC) &F77_NAME(code2), 5},
    {"code3", (DL_FUNC) &F77_NAME(code3), 4},
    {NULL, NULL, 0}
};

void R_init_HDcpDetect(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

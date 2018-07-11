#ifndef MAT_ZZ_P_EXTRA__H
#define MAT_ZZ_P_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/ZZ_p.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* magma output, without assign, using variable var           */
/*------------------------------------------------------------*/
void magma_output(const Mat<ZZ_p>& a, const string & var);

/*------------------------------------------------------------*/
/* magma assignment, using variable var, to name              */
/*------------------------------------------------------------*/
void magma_assign(const Mat<ZZ_p>& a, const string & var, const string & name);

/*------------------------------------------------------------*/
/* random matrix                                              */
/*------------------------------------------------------------*/
void random(Mat<ZZ_p>& a, long n, long m);

#endif

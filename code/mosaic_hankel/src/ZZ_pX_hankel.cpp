#include <NTL/ZZ_pX.h>
#include <NTL/matrix.h>

#include "ZZ_pX_mosaic_hankel.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Hankel matrices                                    */
/* stored as                                          */
/*       a5 a4 a3 a2                                  */
/*       a4 a3 a2 a1                                  */
/*       a3 a2 a1 a0                                  */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void to_dense(Mat<ZZ_p>& Mdense, const ZZ_p_hankel& M){
  long n = M.NumRows();
  long m = M.NumCols();
  Mdense.SetDims(n, m);

  for (long i = 0; i < n; i++)
    for (long j = 0; j < m; j++)
      Mdense[i][j] = M(i,j);
}


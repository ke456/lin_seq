#include <NTL/matrix.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>

#include "mat_ZZ_p_extra.h"
#include "magma_output.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* magma output, without assign, using variable var           */
/*------------------------------------------------------------*/
void magma_output(const Mat<ZZ_p>& a, const string & var){
  cout << "Matrix([";
  for (long i = 0; i < a.NumRows(); i++){
    cout << "[";
    for (long j = 0; j < a.NumCols(); j++){
      cout << a[i][j];
      if (j < a.NumCols()-1)
	cout << ", ";
    }
    cout << "]";
    if (i < a.NumRows()-1)
      cout << ", ";
  }
  cout << "])";
}

/*------------------------------------------------------------*/
/* magma assignment, using variable var, to name              */
/*------------------------------------------------------------*/
void magma_assign(const Mat<ZZ_p>& a, const string & var, const string & name){
  cout << name << ":=";
  magma_output(a, var);
  cout << ";\n";
}

/*------------------------------------------------------------*/
/* random matrix of a given degree                            */
/*------------------------------------------------------------*/
void random(Mat<ZZ_p>& a, long n, long m){
  a.SetDims(n, m);
  for (long i = 0; i < n; i++)
    for (long j = 0; j < m; j++)
      a[i][j] = random_ZZ_p();
}


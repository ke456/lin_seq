#include <NTL/mat_ZZ_p.h>
#include <NTL/ZZ_pX.h>

#include "vec_ZZ_p_extra.h"
#include "ZZ_p_toeplitz.h"
#include "ZZ_p_cauchy_geometric.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Cauchy matrices on geometric progressions          */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/*---------------------------------------------------*/
/* default constructor                               */
/*---------------------------------------------------*/
ZZ_p_cauchy_geometric::ZZ_p_cauchy_geometric() {
  m = 0;
  n = 0;
}

/*---------------------------------------------------*/
/* constructor                                       */
/*---------------------------------------------------*/
ZZ_p_cauchy_geometric::ZZ_p_cauchy_geometric(const ZZ_p& a1, const ZZ_p& b1, const ZZ_p& q, long mm, long nn){
  u1 = a1;
  v1 = b1;
  rho = q;
  sqrt_rho = 0;
  m = mm;
  n = nn;

  if (m == 0 && n == 0){
    t = ZZ_p_toeplitz();
    powers_irho.SetLength(0);
  }
  else{
    prepare_inverses_cauchy(vec_toeplitz, u1, v1, rho, m, n); 
    inverse_powers(powers_irho, rho, m);
    t = ZZ_p_toeplitz(vec_toeplitz, m, n);
  }
}

/*---------------------------------------------------*/
/* dimensions                                        */
/*---------------------------------------------------*/
long ZZ_p_cauchy_geometric::NumRows() const {
  return m;
}

long ZZ_p_cauchy_geometric::NumCols() const {
  return n;
}

/*---------------------------------------------------*/
/* computes output = M*input                         */
/*---------------------------------------------------*/
void ZZ_p_cauchy_geometric::mul_right(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const {
  t.mul_right(output, input);
  long m = NumRows();
  for (long i = 0; i < m; i++)
    output[i] *= powers_irho[i];
}

/*---------------------------------------------------*/
/* computes output = M*input                         */
/*---------------------------------------------------*/
void ZZ_p_cauchy_geometric::mul_right(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const {
  t.mul_right(output, input);
  long m = NumRows();
  long a = input.NumCols();
  for (long i = 0; i < m; i++)
    for (long j = 0; j < a; j++)
      output[i][j] *= powers_irho[i];
}

/*---------------------------------------------------*/
/* computes output = M*input without the diagonal    */
/*---------------------------------------------------*/
void ZZ_p_cauchy_geometric::mul_right_simple(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const {
  t.mul_right(output, input);
}

/*---------------------------------------------------*/
/* computes output = M*input without the diagonal    */
/*---------------------------------------------------*/
void ZZ_p_cauchy_geometric::mul_right_simple(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const {
  t.mul_right(output, input);
}

/*---------------------------------------------------*/
/* computes output = M^t*input                       */
/*---------------------------------------------------*/
void ZZ_p_cauchy_geometric::mul_left(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const {
  Vec<ZZ_p> new_in = input;
  long m = NumRows();
  for (long i = 0; i < m; i++)
    new_in[i] *= powers_irho[i];
  t.mul_left(output, new_in);
}

/*---------------------------------------------------*/
/* computes output = M^t*input                       */
/*---------------------------------------------------*/
void ZZ_p_cauchy_geometric::mul_left(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const {
  Mat<ZZ_p> new_in = input;
  long m = NumRows();
  long a = input.NumCols();
  for (long i = 0; i < m; i++){
    ZZ_p *elts = new_in[i].elts();
    for (long j = 0; j < a; j++)
      elts[j] *= powers_irho[i];
  }
  t.mul_left(output, new_in);
}

/*---------------------------------------------------*/
/* computes output = M^t*input without the diagonal  */
/*---------------------------------------------------*/
void ZZ_p_cauchy_geometric::mul_left_simple(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const {
  t.mul_left(output, input);
}

/*---------------------------------------------------*/
/* computes output = M^t*input without the diagonal  */
/*---------------------------------------------------*/
void ZZ_p_cauchy_geometric::mul_left_simple(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const {
  t.mul_left(output, input);
}

/*---------------------------------------------------*/
/* M as a dense matrix                               */
/*---------------------------------------------------*/
void ZZ_p_cauchy_geometric::to_dense(Mat<ZZ_p>& M) const {
  long m = NumRows();
  long n = NumCols();
  M.SetDims(m, n);
  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++)
      M[i][j] = to_ZZ_p(1) / (u1*power(rho, i)-v1*power(rho, j));
}

/*---------------------------------------------------*/
/* computes                                          */
/* 1/(u1-v1 rho^(-m+1)) ... 1/(u1-v1 rho^(n-1))      */
/* these are the entries of the toeplitz matrix      */
/* (with m rows and n columns)                       */
/*---------------------------------------------------*/
void prepare_inverses_cauchy(Vec<ZZ_p>& inverses, const ZZ_p& u1, const ZZ_p& v1, const ZZ_p& rho, long m, long n){
  
  Vec<ZZ_p> vec_den;
  ZZ_p irho = 1/rho;
  vec_den.SetLength(m+n-1);
  if (m+n-1 == 0)
    return;
  vec_den[0] = -v1*power(irho, m-1);
  for (long i = 1; i < m+n-1; i++){
    vec_den[i] = vec_den[i-1]*rho;
    vec_den[i-1] += u1;
  }
  vec_den[n+m-2] += u1;
  inv(inverses, vec_den);
}

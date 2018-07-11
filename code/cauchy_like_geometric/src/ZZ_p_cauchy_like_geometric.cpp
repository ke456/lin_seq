#include <NTL/mat_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <assert.h>

#include "ZZ_p_toeplitz.h"
#include "ZZ_p_cauchy_geometric.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Cauchy-like matrices on geometric progressions     */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/*---------------------------------------------------*/
/* default constructor                               */
/*---------------------------------------------------*/
ZZ_p_cauchy_like_geometric::ZZ_p_cauchy_like_geometric(){ 
}
  
/*---------------------------------------------------*/
/* constructor                                       */
/*---------------------------------------------------*/
ZZ_p_cauchy_like_geometric::ZZ_p_cauchy_like_geometric(const Mat<ZZ_p>& U, const Mat<ZZ_p>& V, 
							 const ZZ_p& a1, const ZZ_p& b1, const ZZ_p& rho){
  G = U;
  H = V;
  C = ZZ_p_cauchy_geometric(a1, b1, rho, G.NumRows(), H.NumRows());
}

/*---------------------------------------------------*/
/* dimensions                                        */
/*---------------------------------------------------*/
long ZZ_p_cauchy_like_geometric::NumRows() const {
  return G.NumRows();
}

long ZZ_p_cauchy_like_geometric::NumCols() const {
  return H.NumRows();
}

long ZZ_p_cauchy_like_geometric::NumGens() const {
  return G.NumCols();
}

/*---------------------------------------------------*/
/* computes output = M*input                         */
/*---------------------------------------------------*/
void ZZ_p_cauchy_like_geometric::mul_right(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const {

  long m = NumRows();
  long n = NumCols();
  long alpha = NumGens();

  Vec<ZZ_p> new_in, new_out;
  new_in.SetLength(n);

  output.SetLength(m);
  for (long i = 0; i < m; i++)
    output[i] = 0;

  for (long i = 0; i < alpha; i++){
    for (long j = 0; j < n; j++)
      new_in[j] = H[j][i]*input[j];
    C.mul_right(new_out, new_in);
    for (long j = 0; j < m; j++)
      output[j] += G[j][i]*new_out[j];
  }
}

/*---------------------------------------------------*/
/* computes output = M*input                         */
/*---------------------------------------------------*/
void ZZ_p_cauchy_like_geometric::mul_right(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const {

  long m = NumRows();
  long n = NumCols();
  long alpha = NumGens();
  long beta = input.NumCols();
  output.SetDims(m, beta);

  if (m == 0 || n == 0){
    for (long i = 0; i < m; i++){
      ZZ_p * elts = output[i].elts();
      for (long j = 0; j < beta; j++)
	elts[j] = 0;
    }
    return;
  }

  Vec<ZZ_p> new_in, new_out, vec_out;
  new_in.SetLength(n);
  vec_out.SetLength(m);

  const ZZ_p_toeplitz * t = &(C.t);

  for (long i = 0; i < beta; i++){
    for (long j = 0; j < m; j++)
      vec_out[j] = 0;
    for (long k = 0; k < alpha; k++){
      for (long j = 0; j < n; j++)
	new_in[j] = H[j][k]*input[j][i];
      t->mul_right(new_out, new_in);
      for (long j = 0; j < m; j++)
	vec_out[j] += G[j][k]*new_out[j];
    }
    const ZZ_p * diag = C.powers_irho.elts();
    for (long j = 0; j < m; j++)
      output[j][i] = diag[j]*vec_out[j];
  }

}


/*---------------------------------------------------*/
/* computes output = M^t*input                       */
/*---------------------------------------------------*/
void ZZ_p_cauchy_like_geometric::mul_left(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const {

  long m = NumRows();
  long n = NumCols();
  long alpha = NumGens();

  Vec<ZZ_p> new_in, new_out;
  new_in.SetLength(m);

  output.SetLength(n);
  for (long i = 0; i < n; i++)
    output[i] = 0;

  for (long i = 0; i < alpha; i++){
    for (long j = 0; j < m; j++)
      new_in[j] = G[j][i]*input[j];
    C.mul_left(new_out, new_in);
    for (long j = 0; j < n; j++)
      output[j] += H[j][i]*new_out[j];
  }
}


/*---------------------------------------------------*/
/* computes output = M^t*input                       */
/*---------------------------------------------------*/
void ZZ_p_cauchy_like_geometric::mul_left(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const {


  long m = NumRows();
  long n = NumCols();
  long alpha = NumGens();
  long beta = input.NumCols();
  output.SetDims(n, beta);

  if (m == 0 || n == 0){
    for (long i = 0; i < n; i++){
      ZZ_p * elts = output[i].elts();
      for (long j = 0; j < beta; j++)
	elts[j] = 0;
    }
    return;
  }

  Vec<ZZ_p> new_in, new_in_2, new_out, vec_out;
  new_in.SetLength(m);
  new_in_2.SetLength(m);
  vec_out.SetLength(n);

  for (long i = 0; i < beta; i++){
    const ZZ_p * diag = C.powers_irho.elts();
    for (long j = 0; j < m; j++)
      new_in_2[j] = input[j][i] * diag[j];
    for (long j = 0; j < n; j++)
      vec_out[j] = 0;
    for (long k = 0; k < alpha; k++){
      for (long j = 0; j < m; j++)
	new_in[j] = G[j][k]*new_in_2[j];
      C.mul_left_simple(new_out, new_in);
      for (long j = 0; j < n; j++)
	vec_out[j] += H[j][k]*new_out[j];
    }
    for (long j = 0; j < n; j++)
      output[j][i] = vec_out[j];
  }
}

/*---------------------------------------------------*/
/* M as a dense matrix                               */
/*---------------------------------------------------*/
void ZZ_p_cauchy_like_geometric::to_dense(Mat<ZZ_p>& M) const {
  long m = NumRows();
  long n = NumCols();

  M = G*transpose(H);
  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++)
      M[i][j] /= (C.u1*power(C.rho, i)-C.v1*power(C.rho, j));
}

#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_p.h>

#include "ZZ_p_toeplitz.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Toeplitz matrices                                  */
/* stored as                                          */
/*       a2 a3 a4 a5                                  */
/*       a1 a2 a3 a4                                  */
/*       a0 a1 a2 a3                                  */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/*----------------------------------------------------*/
/* default constructor                                */
/*----------------------------------------------------*/
ZZ_p_toeplitz::ZZ_p_toeplitz(){
  data.SetLength(0);
  n = m = 0;
}

/*----------------------------------------------------*/
/* input vector is as showed above                    */
/*----------------------------------------------------*/
ZZ_p_toeplitz::ZZ_p_toeplitz(const Vec<ZZ_p>& input, long rows, long cols){
  n = rows;
  m = cols;
  data.rep.SetLength(n+m-1);
  data_rev.rep.SetLength(n+m-1);
  
  for (long i = 0; i < n+m-1; i++){
    data_rev[i] = input[n+m-2-i];
    data[i] = input[i];
  }
  data.normalize();
  data_rev.normalize();
}

/*---------------------------------------------------*/
/* dimensions                                        */
/*---------------------------------------------------*/
long ZZ_p_toeplitz::NumRows() const {
  return n;
}

long ZZ_p_toeplitz::NumCols() const {
  return m;
}

/*---------------------------------------------------*/
/* data access                                       */
/*---------------------------------------------------*/
const ZZ_p& ZZ_p_toeplitz::operator ()(long i, long j) const {
  return data[n-1-i+j];
}

/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void ZZ_p_toeplitz::to_dense(Mat<ZZ_p>& Mdense) const { 
  long n = NumRows();
  long m = NumCols();
  Mdense.SetDims(n, m);
  
  for (long i = 0; i < n; i++)
    for (long j = 0; j < m; j++)
      Mdense[i][j] = (*this)(i, j);
}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
void ZZ_p_toeplitz::mul_right(Vec<ZZ_p>& res, const Vec<ZZ_p>& input) const {

  long nM = NumRows();
  long mM = NumCols();
  res.SetLength(nM);

  if (nM == 0 || mM == 0){
    for (long i = 0; i < nM; i++)
      res[i] = 0;
    return;
  }

  ZZ_pX input_X;
  input_X.rep = input;
  input_X.normalize();

  ZZ_pX tmp = data_rev * input_X;
  for (long i = 0; i < nM; i++)
    res[i] = coeff(tmp, i + mM - 1);

}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
void ZZ_p_toeplitz::mul_right(Mat<ZZ_p>& res, const Mat<ZZ_p>& input) const {

  long m = NumRows();
  long n = NumCols();

  Vec<ZZ_p> vec_in, vec_out;
  vec_in.SetLength(n);
  long a = input.NumCols();
  res.SetDims(m, a);
  for (long i = 0; i < a; i++){
    ZZ_p *elts = vec_in.elts();
    for (long j = 0; j < n; j++)
      elts[j] = input[j][i];
    mul_right(vec_out, vec_in);
    elts = vec_out.elts();
    for (long j = 0; j < m; j++)
      res[j][i] = elts[j];
  }

}

/*----------------------------------------------------*/
/* left multiplication                                */
/*----------------------------------------------------*/
void ZZ_p_toeplitz::mul_left(Vec<ZZ_p>& res, const Vec<ZZ_p>& input) const {

  long nM = NumRows();
  long mM = NumCols();
  res.SetLength(mM);

  if (nM == 0 || mM == 0){
    for (long i = 0; i < mM; i++)
      res[i] = 0;
    return;
  }

  ZZ_pX input_X;
  input_X.rep = input;
  input_X.normalize();

  ZZ_pX tmp = data * input_X;
  for (long i = 0; i < mM; i++)
    res[i] = coeff(tmp, i + nM - 1);
}

/*----------------------------------------------------*/
/* left multiplication                                */
/*----------------------------------------------------*/
void ZZ_p_toeplitz::mul_left(Mat<ZZ_p>& res, const Mat<ZZ_p>& input) const {
  long m = NumRows();
  long n = NumCols();

  Vec<ZZ_p> vec_in, vec_out;
  vec_in.SetLength(m);
  long a = input.NumCols();
  res.SetDims(n, a);
  for (long i = 0; i < a; i++){
    ZZ_p *elts = vec_in.elts();
    for (long j = 0; j < m; j++)
      elts[j] = input[j][i];
    mul_left(vec_out, vec_in);
    elts = vec_out.elts();
    for (long j = 0; j < n; j++)
      res[j][i] = elts[j];
  }
}

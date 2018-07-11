#include <NTL/tools.h>
#include <time.h>
#include <cmath>

#include "ZZ_p_block_sylvester.h"
#include "ZZ_pX_extra.h"
#include "mat_ZZ_pX_extra.h"

NTL_CLIENT

/*----------------------------------------------------*/
/* is this useful?                                    */
/*----------------------------------------------------*/
void ZZ_p_bivariate_modular_composition::init (const ZZ_pX &f, const Vec<long> &type_new, long prec_new){
  initialized = true;
  type = type_new;
  f_field = f;
  
  prec = prec_new;
  sqrtP = ceil(sqrt(type.length()));
  ZZ_pX running = ZZ_pX{1};
  Vec<ZZ_pX> fs;
  
  fs.SetLength(sqrtP);
  Vec<long> Stype;
  Stype.SetLength(sqrtP);

  ZZ_pX_poly_multiplier multF(f, prec);
  for (long i = 0; i < sqrtP; i++){
    Stype[i] = max_of_type;
    fs[i] = running;
    multF.mul(running, running);
    trunc(running, running, prec);
  }

  F_field = running;

  S = ZZ_p_block_sylvester_general(fs, Stype, prec);
}

/*----------------------------------------------------*/
/* multiplies rhs using Horner's rule                 */
/*----------------------------------------------------*/
Vec<ZZ_p> ZZ_p_bivariate_modular_composition::mul_right_Horners(const Vec<ZZ_p> &rhs){
  if (!initialized) 
    throw "must init first";
  Vec<ZZ_pX> rhs_poly;
  create_lhs_list(rhs_poly, rhs);
  ZZ_pX result = rhs_poly[rhs_poly.length()-1];
  for (long i = rhs_poly.length()-2; i >= 0; i--)
    result = trunc((result * f_field), prec) + rhs_poly[i];

  Vec<ZZ_p> v;
  v.SetLength(prec);
  for (long i = 0; i < prec; i++)
    v[i] = coeff(result, i);
  return v;
}

/*----------------------------------------------------*/
/* mult using the baby steps / giant steps algorithm  */
/*----------------------------------------------------*/
Vec<ZZ_p> ZZ_p_bivariate_modular_composition::mul_left(const Vec<ZZ_p> &lhs){

  if (!initialized) 
    throw "must init first";

  Mat<ZZ_p> lhs_mat;
  lhs_mat.SetDims(prec, ceil((type.length()*1.0) / sqrtP));
  
  ZZ_pX tmp;
  tmp.rep.SetLength(prec);
  for (long j = 0; j < prec; j++){
    lhs_mat[j][0] = lhs[j];
    tmp.rep[j] = lhs[prec - 1 - j];
  }
  tmp.normalize();

  for (long i = 1; i < lhs_mat.NumCols(); i++){
    tmp = trunc(F_field * tmp, prec);
    for (long j = 0; j < prec; j++)
      lhs_mat[j][i] = coeff(tmp, prec - 1 - j);
  }

  Mat<ZZ_p> prod = S.mul_left(lhs_mat);

  Vec<ZZ_p> result;
  result.SetLength(num_cols);
  long idx = 0;
  long nb = 0;
  for (long u = 0; u < prod.NumCols(); u++){
    long acc = 0;
    for (long v = 0; v < sqrtP; v++)
      if (idx < type.length()){
	for (long j = 0; j < type[idx] + 1; j++)
	  if (nb < num_cols){
	    result[nb] = prod[acc + j][u];
	    nb++;
	  }
	idx++;
	acc += max_of_type + 1;
      }
  }

  return result;
}

/*----------------------------------------------------*/
/* mult using the baby steps / giant steps algorithm  */
/*----------------------------------------------------*/
Vec<ZZ_p> ZZ_p_bivariate_modular_composition::mul_right_comp(const Vec<ZZ_p> &rhs){

  if (!initialized) 
    throw "must init first";

  Mat<ZZ_p> rhs_mat;
  rhs_mat.SetDims(sqrtP * (max_of_type+1), ceil((type.length()*1.0) / sqrtP));
  long idx = 0;
  long nb = 0;
  for (long k = 0; k < rhs_mat.NumCols(); k++){
    long acc = 0;
    for (long i = 0; i < sqrtP; i++){
      if (idx < type.length()){
	for (long j = 0; j < type[idx]+1; j++)
	  if (nb < num_cols){
	    rhs_mat[acc + j][k] = rhs[nb];
	    nb++;
	  }
	acc += max_of_type + 1;
	idx++;
      }
    }
  }

  Mat<ZZ_p> prod = S.mul_right(rhs_mat);

  Vec<ZZ_pX> B1;
  B1.SetLength(prod.NumCols());
  for (long i = 0; i < prod.NumCols(); i++){
    ZZ_pX a;
    a.rep.SetLength(prec);
    for (long j = 0; j < prec; j++)
      a.rep[j] = prod[j][i];
    a.normalize();
    B1[i] = a;
  }

  ZZ_pX p;
  p = B1[B1.length() - 1];
  ZZ_pX_poly_multiplier multF(F_field, prec);
  for (long i = B1.length()-2; i >= 0; i--){
    multF.mul(p, p);
    p = trunc(p, prec) + B1[i];
  }

  Vec<ZZ_p> v;
  v.SetLength(prec);
  for (long i = 0; i < prec; i++)
    v[i] = coeff(p, i);

  return v;
}

/*----------------------------------------------------*/
/* header for multiplying                             */
/* TODO: thresholds                                   */
/*----------------------------------------------------*/
Vec<ZZ_p> ZZ_p_bivariate_modular_composition::mul_right(const Vec<ZZ_p> &rhs){
  return mul_right_comp(rhs);
}

/*----------------------------------------------------*/
/* makes a dense matrix                               */
/*----------------------------------------------------*/
void ZZ_p_bivariate_modular_composition::to_dense(Mat<ZZ_p> & dense){
  dense.SetDims(prec, num_cols);

  long idx = 0;
  ZZ_pX tmp;
  SetCoeff(tmp, 0, to_ZZ_p(1));
  for (long i = 0; i < type.length(); i++){
    for (long k = 0; k < type[i] + 1; k++){
      for (long j = k; j < prec; j++)
	dense[j][idx] = coeff(tmp, j - k);
      idx++;
    }
    tmp = trunc(tmp * f_field, prec);
  }
}


/*----------------------------------------------------*/
/* default constructor; does nothing                  */
/*----------------------------------------------------*/
ZZ_p_bivariate_modular_composition::ZZ_p_bivariate_modular_composition(){
}

/*----------------------------------------------------*/
/* input: Vec of polynomials fs                       */
/*        type                                        */
/*        output precision                            */
/*----------------------------------------------------*/
ZZ_p_bivariate_modular_composition::ZZ_p_bivariate_modular_composition(const ZZ_pX& f, 
								       const Vec<long> &type, long prec): 
  ZZ_p_block_sylvester(type, prec) {
  init(f, type, prec);
}

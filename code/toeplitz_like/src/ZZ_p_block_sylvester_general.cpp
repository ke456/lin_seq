#include "ZZ_pX_extra.h"
#include "mat_ZZ_pX_extra.h"
#include "ZZ_p_block_sylvester.h"

NTL_CLIENT

/*----------------------------------------------------*/
/* sets up the matrices                               */
/*----------------------------------------------------*/
void ZZ_p_block_sylvester_general::init_mat(){

  // this matrix holds chunks of the fs for right-mul
  long rows = ceil( (1.0*prec) / (1.0*(max_of_type+1)));
  matF_ZZ.SetDims(rows, type.length());
  
  for (long j = 0; j < matF_ZZ.NumCols(); j++){
    long acc = 0;
    for (long i = 0; i < matF_ZZ.NumRows(); i++){
      ZZX p;
      for (long s = 0; s < max_of_type + 1; s++)
	SetCoeff(p, s, coeff(f_ZZ[j], acc + s));
	// SetCoeff(p, s, conv<ZZ>(coeff(fs[j], acc + s)));
      matF_ZZ.put(i, j, p);
      acc += max_of_type + 1;
    }
  }

  // this matrix holds chunks of the fs for left-mul
  matF_left_ZZ.SetDims(rows, type.length());

  for (long j = 0; j < matF_ZZ.NumCols(); j++){
    long acc = 0;
    for (long i = 0; i < matF_ZZ.NumRows(); i++){
      ZZX p;
      for (long s = 0; s < 2*(max_of_type + 1)-1; s++)
	if (acc + max_of_type - s >= 0)
	  SetCoeff(p, s, coeff(f_ZZ[j], acc + max_of_type - s));
	  // SetCoeff(p, s, conv<ZZ>(coeff(fs[j], acc + max_of_type - s)));
      matF_left_ZZ.put(i, j, p);
      acc += max_of_type + 1;
    }
  }
  mat_init = true;
}


/*----------------------------------------------------*/
/*----------------------------------------------------*/
void ZZ_p_block_sylvester_general::init(const Vec<ZZ_pX>& fs, const Vec<long> &type, long prec){
  this->type = type;
  this->prec = prec;
  f_ZZ = conv<Vec<ZZX>>(fs);
  
  f_rev_ZZ.SetLength(type.length());
  for (long i = 0; i < type.length(); i++){
    f_rev_ZZ[i].rep.SetLength(prec);
    for (long j = 0; j < prec; j++)
      f_rev_ZZ[i].rep[j] = conv<ZZ>(coeff(fs[i], prec-1-j));
    f_rev_ZZ[i].normalize();
  }

  initialized = true;
}


/*----------------------------------------------------*/
/*----------------------------------------------------*/
void ZZ_p_block_sylvester_general::init(const Vec<ZZX>& fs, const Vec<long> &type, long prec){
  this->type = type;
  this->prec = prec;
  f_ZZ = fs;
  
  f_rev_ZZ.SetLength(type.length());
  for (long i = 0; i < type.length(); i++){
    f_rev_ZZ[i].rep.SetLength(prec);
    for (long j = 0; j < prec; j++)
      f_rev_ZZ[i].rep[j] = coeff(fs[i], prec-1-j);
    f_rev_ZZ[i].normalize();
  }

  initialized = true;
}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
Vec<ZZ_p> ZZ_p_block_sylvester_general::mul_right(const Vec<ZZ_p> &rhs){
  // maybe better?
  if (!initialized) 
    throw "must call init first";
  Vec<ZZ_pX> rhs_poly;
  create_lhs_list(rhs_poly, rhs);
  Vec<ZZX> rhs_ZZ = conv<Vec<ZZX>>(rhs_poly);
  ZZX result_ZZ;
  for (long i = 0; i < f_ZZ.length(); i++)
    result_ZZ += trunc(f_ZZ[i] * rhs_ZZ[i], prec);
  Vec<ZZ_p> v;
  v.SetLength(prec);
  for (long i = 0; i < prec; i++)
    v[i] = conv<ZZ_p>(coeff(result_ZZ, i));
  return v;
  // Vec<ZZ_pX> f = conv<Vec<ZZ_pX>>(f_ZZ);
  // if (!initialized) 
  //   throw "must call init first";
  // Vec<ZZ_pX> rhs_poly;
  // create_lhs_list(rhs_poly, rhs);
  // ZZ_pX result;
  // for (long i = 0; i < f.length(); i++)
  //   result += trunc(f[i] * rhs_poly[i], prec);
  // Vec<ZZ_p> v;
  // v.SetLength(prec);
  // for (long i = 0; i < prec; i++)
  //   v[i] = coeff(result, i);
  // return v;
}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
Mat<ZZ_p> ZZ_p_block_sylvester_general::mul_right(const Mat<ZZ_p> &rhs){

  if (!mat_init)
    init_mat();

  Mat<ZZ_pX> matF = conv<Mat<ZZ_pX>>(matF_ZZ);
  Mat<ZZ_pX> rhs_poly_mat;
  rhs_poly_mat.SetDims(type.length(), rhs.NumCols());
  Vec<ZZ_p> column;
  column.SetLength(rhs.NumRows());

  for (long i = 0; i < rhs.NumCols(); i++){
    for (long j = 0; j < rhs.NumRows(); j++)
      column[j] = rhs[j][i];
    Vec<ZZ_pX> rhs_poly;
    create_lhs_list(rhs_poly, column);
    for (long j = 0; j < type.length(); j++)
      rhs_poly_mat[j][i] = rhs_poly[j];
  }

  Mat<ZZ_pX> product;
  mul_CRT_CTFT(product, matF, rhs_poly_mat);

  Mat<ZZ_p> result;
  result.SetDims(prec, rhs.NumCols());
  
  for (long  i = 0; i < result.NumCols(); i++){

    const ZZ_p * cC = product[0][i].rep.elts();
    long old_len = product[0][i].rep.length();
    for (long j = 0; j < min(old_len, prec); j++)
      result[j][i] = cC[j];
    long idx = max_of_type + 1;

    for (long a = 1; a < product.NumRows(); a++){
      const ZZ_p * cC = product[a][i].rep.elts();
      long new_len = product[a][i].rep.length();
      long b1 = min(min(old_len-(max_of_type+1), new_len), prec-idx);
      long b2 = min(new_len, prec-idx);
      for (long j = 0; j < b1; j++)
	      result[j+idx][i] += cC[j];
      for (long j = max(0, old_len-(max_of_type+1)); j < b2; j++)
	      result[j+idx][i] = cC[j];
      old_len = new_len;
      idx += max_of_type + 1;
    }
  }

  return result;
}

/*----------------------------------------------------*/
/* left multiplication                                */
/*----------------------------------------------------*/
Vec<ZZ_p> ZZ_p_block_sylvester_general::mul_left(const Vec<ZZ_p> &in){
  if (!initialized) 
    throw "must call init first";

  ZZ_pX in_poly;
  in_poly.rep = in;
  in_poly.normalize();
  // Vec<ZZ_pX> f = conv<Vec<ZZ_pX>>(f_ZZ);
  Vec<ZZ_pX> f_rev = conv<Vec<ZZ_pX>>(f_rev_ZZ);

  Vec<ZZ_p> out;
  out.SetLength(num_cols);
  long acc = 0;
  for (long i = 0; i < f_ZZ.length(); i++){
    ZZ_pX result = in_poly * f_rev[i];
    for (long j = 0; j < type[i] + 1; j++)
      out[acc++] = coeff(result, prec - 1 + j);
  }
  return out;
}

/*----------------------------------------------------*/
/* left multiplication                                */
/*----------------------------------------------------*/
Mat<ZZ_p> ZZ_p_block_sylvester_general::mul_left(const Mat<ZZ_p> &in){

  
  if (!mat_init)
    init_mat();

  Mat<ZZ_pX> matF_left = conv<Mat<ZZ_pX>>(matF_left_ZZ);

  Mat<ZZ_pX> poly_in;
  poly_in.SetDims(in.NumCols(), matF_left.NumRows());
  for (long i = 0; i < in.NumCols(); i++){
    long acc = 0;
    for (long k = 0; k < matF_left.NumRows(); k++){
      ZZ_pX tmp;
      for (long j = 0; j < max_of_type + 1; j++)
	if (acc < in.NumRows())
	  SetCoeff(tmp, j, in[acc++][i]);
      poly_in[i][k] = tmp;
    }
  }

  Mat<ZZ_pX> product;
  mul_CRT_CTFT(product, poly_in, matF_left);

  Mat<ZZ_p> res;
  res.SetDims(num_cols, in.NumCols());
  for (long i = 0; i < in.NumCols(); i++){
    long acc = 0;
    for (long j = 0; j < product.NumCols(); j++)
      for (long k = 0; k < max_of_type + 1; k++)
	res[acc++][i] = coeff(product[i][j], k + max_of_type);
  }

  return res;

}

/*----------------------------------------------------*/
/* makes a dense matrix                               */
/*----------------------------------------------------*/
void ZZ_p_block_sylvester_general::to_dense(Mat<ZZ_p> & dense){
  dense.SetDims(prec, num_cols);
  long acc = 0;
  Vec<ZZ_pX> f = conv<Vec<ZZ_pX>>(f_ZZ);
  for (long i = 0; i < type.length(); i++){
    for (long k = 0; k < type[i]+1; k++){
      for (long ell = 0; ell < prec - k; ell++)
	dense[ell + k][acc] = coeff(f[i], ell);
      acc++;
    }
  }

};


/*----------------------------------------------------*/
/* input: Vec of polynomials fs                       */
/*        type                                        */
/*        output precision                            */
/*----------------------------------------------------*/
ZZ_p_block_sylvester_general::ZZ_p_block_sylvester_general(const Vec<ZZ_pX>& fs, const Vec<long> &type, long prec):
  ZZ_p_block_sylvester(type, prec) {
  init(fs, type, prec);
}


/*----------------------------------------------------*/
/* input: Vec of polynomials fs                       */
/*        type                                        */
/*        output precision                            */
/*----------------------------------------------------*/
ZZ_p_block_sylvester_general::ZZ_p_block_sylvester_general(const Vec<ZZX>& fs, const Vec<long> &type, long prec):
  ZZ_p_block_sylvester(type, prec) {
  init(fs, type, prec);
}

/*----------------------------------------------------*/
/* does nothing                                       */  
/*----------------------------------------------------*/
ZZ_p_block_sylvester_general::ZZ_p_block_sylvester_general():
  ZZ_p_block_sylvester() {

}

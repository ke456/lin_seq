#include <cstdlib>
#include <cmath>

#include "magma_output.h"
#include "ZZ_hermite_pade.h"
#include "vec_ZZ_p_extra.h"
#include "lzz_pXY.h"
#include "lzz_p_extra.h"

NTL_CLIENT

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/
/* hermite pade for algebraic approximants                        */
/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/

/*----------------------------------------------------------------*/
/* creates a new bivariate modular composition                    */
/*----------------------------------------------------------------*/
SmartPtr<ZZ_p_block_sylvester> hermite_pade_algebraic::create_bmc(){
  ZZ_pX f_p;
  ZZ_p denom_p= conv<ZZ_p>(denom);
  ZZ_p i_denom = 1 / denom_p;;
  conv(f_p, f_full_prec);
  f_p /= denom_p;
  return MakeSmart<ZZ_p_bivariate_modular_composition>(ZZ_p_bivariate_modular_composition(f_p, type, sizeX));
}

/*----------------------------------------------------------------*/
/* gives a candidate for a reduced type by GCD techniques         */
/*----------------------------------------------------------------*/
Vec<long> hermite_pade_algebraic::find_new_type(){
  if (rank + 1 == sizeY) 
    return type;

  switch_context(0);
	
  Vec<zz_pX> v1, v2;
  hermite_pade::random_solution_mod_p(v1);
  hermite_pade::random_solution_mod_p(v2);
	
  // magma_init();
  // magma_assign(v1, "v1");
  // magma_assign(v2, "v2");
  zz_pXY b1(v1);
  zz_pXY b2(v2);
  zz_pXY gcd;
  GCD(gcd, b1, b2);

  // extracting the type
  Vec<long> t;
  for (long i = 0; i < gcd.coeffX.length(); i++)
    t.append(deg(gcd.coeffX[i]));
  return t;
}


/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/
void hermite_pade_algebraic::generators_hankel_last_row_last_column(Mat<ZZ_p> & X, Mat<ZZ_p> & Y, Vec<ZZ_p> & r, Vec<ZZ_p> & c){
  
  // setting up the mosaic Hankel matrix
  Vec<ZZ_p_hankel> vec_H;
  ZZ_pX running {1};
  ZZ_p denom_p = conv<ZZ_p>(denom);
  ZZ_pX f = conv<ZZ_pX>(f_full_prec);

  for (long i = 0; i < type.length(); i++){
    Vec<ZZ_p> v;
    Vec<ZZ_p> inp_vec;
    conv(v, running);
    v.SetLength(prec, ZZ_p{0});
    for (int j = 0; j < v.length(); j++)
      inp_vec.append(v[v.length() - 1 - j]);
    for (int j = 0; j < type[i]; j++)
      inp_vec.append(ZZ_p{0});
    vec_H.append(ZZ_p_hankel(inp_vec, prec, type[i] + 1));
    running = trunc(running * f, prec);
    running /= denom_p;
  }
  Vec<Vec<ZZ_p_hankel>> hankel_matrices;
  hankel_matrices.append(vec_H);
  ZZ_p_mosaic_hankel MH(hankel_matrices);
  generators(X, Y, MH);
  last_row_of_block(r, MH.NumBlockRows()-1, MH);
  last_column_of_block(c, MH.NumBlockCols()-1, MH);
}



/*----------------------------------------------------------------*/
/* initializes everything                                         */
/*----------------------------------------------------------------*/
void hermite_pade_algebraic::init(const Vec<long> &type_in){
  level = 0;
  ZZ_p::init(p_powers[0]);
  zz_pX f_field = conv<zz_pX>(f_full_prec);
  type = type_in;

  // setting up the mosaic Hankel matrix
  Vec<hankel> vec_H;
  zz_pX running {1};
  zz_p denom_p = conv<zz_p>(denom);

  for (long i = 0; i < type.length(); i++){
    Vec<zz_p> v;
    Vec<zz_p> inp_vec;
    conv(v, running);
    v.SetLength(prec, zz_p{0});
    for (int j = 0; j < v.length(); j++)
      inp_vec.append(v[v.length() - 1 - j]);
    for (int j = 0; j < type[i]; j++)
      inp_vec.append(zz_p{0});
    vec_H.append(hankel(inp_vec, prec, type[i] + 1));
    running = trunc(running * f_field, prec);
    running /= denom_p;
  }
  Vec<Vec<hankel>> hankel_matrices;
  hankel_matrices.append(vec_H);
  mosaic_hankel MH(hankel_matrices);

  lzz_p_cauchy_like_geometric CL; // the cauchy matrix
  zz_pX_Multipoint_Geometric X_int, Y_int; // preconditioners
  
  // find the rank of the original matrix
  Vec<zz_p> e_zz_p, f_zz_p;                            // the diagonal matrices
  to_cauchy_grp(CL, X_int, Y_int, e_zz_p, f_zz_p, MH); // converting from Hankel to Cauchy

  rank = invert(invA, CL);                             // inverting M mod p
  // if (rank == -1){
  //   cout << type << endl;
  //   Mat<zz_p> CLD;
  //   CL.to_dense(CLD);
  //   cout << "k:=GF(" << zz_p::modulus() << ";\n";
  //   magma_assign(CLD, "cld");
  // }

  sizeX = X_int.length();
  sizeY = Y_int.length();
  	
  // reset all vectors
  vec_M.SetLength(0);
  vec_X_int.SetLength(0);
  vec_Y_int.SetLength(0);

  // sets up the first M
  set_up_bmc();

  // converting the preconditioners that do not change
  this->e = conv<Vec<ZZ>>(e_zz_p);
  this->f = conv<Vec<ZZ>>(f_zz_p);
  zz_p c_zz, d_zz; 
  X_int.point(c_zz, 0);
  Y_int.point(d_zz, 0); 
  c = conv<ZZ>(c_zz);
  d = conv<ZZ>(d_zz);

  // initializing the X_int and Y_int stuff
  zz_p w_zz_p, w2;
  X_int.point(w_zz_p, 1);
  Y_int.point(w2, 1);
  w_zz_p = w_zz_p / c_zz;
  this->w = w_zz_p.LoopHole();
  ZZ_p w_p = conv<ZZ_p>(this->w);
  // find the order of w
  order = order_dyadic(w_zz_p);

  vec_X_int.append(ZZ_pX_Multipoint_FFT(w_p, conv<ZZ_p>(this->c), sizeX));
  vec_Y_int.append(ZZ_pX_Multipoint_FFT(w_p, conv<ZZ_p>(this->d), sizeY));
}


/*----------------------------------------------------------------*/
/* f: the power series                                            */
/* type: the type of the approximant                              */
/* sigma: requested precision                                     */
/*----------------------------------------------------------------*/
hermite_pade_algebraic::hermite_pade_algebraic(const ZZX &f, const Vec<long>& type, long sigma, long fft_index):
  hermite_pade(fft_index) {
  stop_criterion = 1;
  f_full_prec = f;
  prec = sigma;
  init(type);
}

/*----------------------------------------------------------------*/
/* f: the power series                                            */
/* type: the type of the approximant                              */
/* sigma: requested precision                                     */
/* w: solution modulo the next fft prime                          */
/*----------------------------------------------------------------*/
hermite_pade_algebraic::hermite_pade_algebraic(const ZZX &f, const Vec<long>& type, long sigma, Vec<long>& w, long fft_index):
  hermite_pade(fft_index) {
  f_full_prec = f;
  prec = sigma;
  init(type);
  witness = w;
  stop_criterion = 2;
}

/*----------------------------------------------------------------*/
/* f: numerator of the power series                               */
/* denom: its denominator                                         */
/* type: the type of the approximant                              */
/* sigma: requested precision                                     */
/*----------------------------------------------------------------*/
hermite_pade_algebraic::hermite_pade_algebraic(const ZZX &f, const ZZ &denom_in, const Vec<long>& type, long sigma, long fft_index):
  hermite_pade(fft_index) {
  stop_criterion = 1;
  f_full_prec = f;
  prec = sigma;
  denom = denom_in;
  init(type);
}

/*----------------------------------------------------------------*/
/* does nothing                                                   */
/*----------------------------------------------------------------*/
hermite_pade_algebraic::hermite_pade_algebraic(){
}

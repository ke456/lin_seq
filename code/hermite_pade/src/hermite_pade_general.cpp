#include "vec_ZZ_p_extra.h"
#include "lzz_pXY.h"
#include "lzz_p_extra.h"
#include "ZZ_hermite_pade.h"
#include "ZZ_p_block_sylvester.h"

NTL_CLIENT

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/
/* hermite pade for general polynomials                           */
/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/

/*---------------------------------------------------------*/
/* tries to increase the rank of M by adding random blocks */
/*---------------------------------------------------------*/
void hermite_pade_general::increase_rank(Vec<hankel> & vh_zz_p, Vec<ZZ_hankel> & vh_ZZ, long add){
  // cout << "adding " << add << endl;
  rows_added = add;
  vh_zz_p.SetLength(0);
  vh_ZZ.SetLength(0);
  
  for (long i = 0; i < type.length(); i++){
    Vec<long> adding;
    Vec<long> values;
    long to_add = max(1, add - type[i]);
    adding.SetLength(add+type[i],0);
    for (long j = 0; j < to_add; j++){
      values.append(rand() % 100);
      adding[j+add-to_add] = values[j];
    }
    vec_added.append(values);

    Vec<zz_p> add_p = conv<Vec<zz_p>>(adding);
    Vec<ZZ> add_ZZ = conv<Vec<ZZ>>(adding);
    hankel h_zz_p{add_p, add, type[i]+1};
    ZZ_hankel h_ZZ{add_ZZ, add, type[i]+1};
    vh_zz_p.append(h_zz_p);
    vh_ZZ.append(h_ZZ);
  }
}

/*----------------------------------------------------------------*/
/* creates a new block Sylvester matrix                           */
/*----------------------------------------------------------------*/
SmartPtr<ZZ_p_block_sylvester> hermite_pade_general::create_bmc(){
  return MakeSmart<ZZ_p_block_sylvester_general>(ZZ_p_block_sylvester_general(vec_fs, type, original_sizeX));
  // return MakeSmart<ZZ_p_block_sylvester_general>(ZZ_p_block_sylvester_general(fs_p, type, original_sizeX));
}

/*----------------------------------------------------------------*/
/* does the product by M and the extra rows                       */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade_general::mul_M_right(const Vec<ZZ_p> &b){
  Vec<ZZ_p> upper = hermite_pade::mul_M_right(b);
  if (rows_added != 0){
    Vec<Vec<ZZ_p>> b_split = split_on_type(flip_on_type(b));
    double time = GetTime();
    ZZ_pX lower;
    for (long i = 0; i < type.length(); i++){
      Vec<ZZ_p> values;
      for (long j = vec_added[i].length()-1; j >= 0; j--)
	values.append(ZZ_p(vec_added[i][j]));
      lower += conv<ZZ_pX>(values) * conv<ZZ_pX>(b_split[i]);
    }
    upper.SetLength(original_sizeX);
    for (long i = 0; i < deg(lower)+1; i++){
      upper.append(lower[i]);
    }
    time_spent_lower += GetTime() - time;
  }
  return upper;
}

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/
void hermite_pade_general::generators_hankel_last_row_last_column(Mat<ZZ_p> & X, Mat<ZZ_p> & Y, Vec<ZZ_p> & r, Vec<ZZ_p> & c){
  X = conv<Mat<ZZ_p>>(G);
  Y = conv<Mat<ZZ_p>>(H);

  Vec<ZZ> row_ZZ;
  last_row_of_block(row_ZZ, MH_ZZ.NumBlockRows()-1, MH_ZZ);
  r = conv<Vec<ZZ_p>>(row_ZZ);

  Vec<ZZ> col_ZZ;
  last_column_of_block(col_ZZ, MH_ZZ.NumBlockCols()-1, MH_ZZ);
  c = conv<Vec<ZZ_p>>(col_ZZ);
}


/*----------------------------------------------------------------*/
/* initializes everything                                         */
/*----------------------------------------------------------------*/
void hermite_pade_general::init(){

  // setting up the mosaic Hankel Matrix over zz_p and ZZ
  Vec<hankel> vec_H;
  Vec<ZZ_hankel> vec_H_ZZ;
  for (long i = 0; i < vec_fs.length(); i++){
    Vec<ZZ> v_ZZ = conv<Vec<ZZ>>(vec_fs[i]);
    Vec<zz_p> v = conv<Vec<zz_p>>(v_ZZ);
    v.SetLength(prec);
    v_ZZ.SetLength(prec);
    Vec<zz_p> inp_zz_p;
    Vec<ZZ> inp_ZZ;
    for (long j = 0; j < v.length(); j++){
      inp_zz_p.append(v[v.length() - 1 - j]);
      inp_ZZ.append(v_ZZ[v_ZZ.length() - 1 - j]);
    }
    for (long j = 0; j < type[i]; j++){
      inp_zz_p.append(zz_p{0});
      inp_ZZ.append(ZZ{0});
    }
    vec_H.append(hankel(inp_zz_p, prec, type[i]+1)); 
    vec_H_ZZ.append(ZZ_hankel(inp_ZZ, prec, type[i]+1)); 
  }

  Vec<Vec<hankel>> hankel_matrices_zz_p;
  hankel_matrices_zz_p.append(vec_H);
  mosaic_hankel MH(hankel_matrices_zz_p);

  Vec<Vec<ZZ_hankel>> hankel_matrices_ZZ;
  hankel_matrices_ZZ.append(vec_H_ZZ);
  
  lzz_p_cauchy_like_geometric CL;                      // the cauchy matrix
  zz_pX_Multipoint_Geometric X_int, Y_int;             // preconditioners
  Vec<zz_p> e_zz_p, f_zz_p;                            // the diagonal matrices

  to_cauchy_grp(CL, X_int, Y_int, e_zz_p, f_zz_p, MH); // converting from Hankel to Cauchy
  rank = invert_fast(invA, CL); // inverting M mod p
  cout << "original rank: " << rank << endl;

  sizeX = X_int.length();
  sizeY = Y_int.length();
  original_sizeX = sizeX;

  if (sizeY -rank != 1){
    Vec<hankel> new_zz_p;
    Vec<ZZ_hankel> new_ZZ;
    increase_rank(new_zz_p, new_ZZ, sizeY-rank-1);
    hankel_matrices_zz_p.append(new_zz_p);
    hankel_matrices_ZZ.append(new_ZZ);
    MH = mosaic_hankel(hankel_matrices_zz_p);
  }
  else{
    rows_added = 0;
  }
  MH_ZZ = ZZ_mosaic_hankel(hankel_matrices_ZZ);
  generators(G, H, MH_ZZ);

  to_cauchy_grp(CL, X_int, Y_int, e_zz_p, f_zz_p, MH); // converting from Hankel to Cauchy
  rank = invert_fast(invA, CL); // inverting M mod p
  cout << "new rank: " << rank << endl;
  sizeX = X_int.length();
  sizeY = Y_int.length();

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

  ZZ_pX_Multipoint_FFT X_int_ZZ_p(w_p, conv<ZZ_p>(this->c), sizeX);
  ZZ_pX_Multipoint_FFT Y_int_ZZ_p(w_p, conv<ZZ_p>(this->d), sizeY);

  vec_w.SetLength(1);
  vec_M.SetLength(0);
  vec_X_int.SetLength(1);
  vec_Y_int.SetLength(1);

  vec_w[0] = ZZ(w);
  set_up_bmc();
  vec_X_int[0] = X_int_ZZ_p;
  vec_Y_int[0] = Y_int_ZZ_p;

}

/*----------------------------------------------------------------*/
/* fs: the power series                                           */
/* type: the type of the approximant                              */
/* sigma: requested precision                                     */
/*----------------------------------------------------------------*/
hermite_pade_general::hermite_pade_general(const Vec<ZZX> &fs, const Vec<long> &type, long sigma, long fft_init): 
  hermite_pade(fft_init){

  prec = sigma;
  this->type = type;
  vec_fs = fs;

  init();
}

/*----------------------------------------------------------------*/
/* fs: the power series                                           */
/* type: the type of the approximant                              */
/* sigma: requested precision                                     */
/* w: solution modulo the next fft prime                          */
/*----------------------------------------------------------------*/
hermite_pade_general::hermite_pade_general(const Vec<ZZX> &fs, const Vec<long> &type, long sigma, Vec<long>& w, long fft_init): 
  hermite_pade(fft_init){
  prec = sigma;
  this->type = type;
  vec_fs = fs;
  init();
  witness = w;
  stop_criterion = 2;
}

/*----------------------------------------------------------------*/
/* does nothing                                                   */
/*----------------------------------------------------------------*/
hermite_pade_general::hermite_pade_general(){
}

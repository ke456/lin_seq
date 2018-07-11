
#include "ZZ_pX_extra.h"
#include "magma_output.h"
#include "ZZ_hermite_pade.h"
#include "vec_ZZ_p_extra.h"
#include "lzz_pXY.h"
#include "lzz_p_extra.h"

NTL_CLIENT

/*----------------------------------------------------------------*/
/* applies a block reversal to v                                  */
/* e.g., type = [1,2] v = [v1,v2,v3,v4,v5]                        */
/* returns [v2,v1,v5,v4,v3] (blocks have length type[i]+1)        */                                    
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::flip_on_type (const Vec<ZZ_p> &v){
  Vec<ZZ_p> r;
  r.SetMaxLength(v.length());
  long acc = 0;
  for (long i = 0; i < type.length(); i++){
    for(long j = type[i]; j >= 0; j--)
      r.append(v[j+acc]);
    acc += type[i] + 1;
  }
  return r;
}


/*----------------------------------------------------------------*/
/* applies a block reversal to v                                  */
/* e.g., type = [1,2] v = [v1,v2,v3,v4,v5]                        */
/* returns [v2,v1,v5,v4,v3] (blocks have length type[i]+1)        */                                    
/*----------------------------------------------------------------*/
Vec<long> hermite_pade::flip_on_type (const Vec<long> &v){
  Vec<long> r;
  r.SetMaxLength(v.length());
  long acc = 0;
  for (long i = 0; i < type.length(); i++){
    for(long j = type[i]; j >= 0; j--)
      r.append(v[j+acc]);
    acc += type[i] + 1;
  }
  return r;
}

/*----------------------------------------------------------------*/
/* applies a block reversal to v                                  */
/*----------------------------------------------------------------*/
Vec<Vec<ZZ>> hermite_pade::flip_on_type (const Vec<Vec<ZZ>> &v){
  Vec<Vec<ZZ>> r;
  r.SetMaxLength(v.length());
  long acc = 0;
  for (long i = 0; i < type.length(); i++){
    for(long j = type[i]; j >= 0; j--)
      r.append(v[j+acc]);
    acc += type[i] + 1;
  }
  return r;
}

/*----------------------------------------------------------------*/
/* splits v according to the current type                         */
/*----------------------------------------------------------------*/
Vec<Vec<ZZ_p>> hermite_pade::split_on_type(const Vec<ZZ_p> &v){
  long acc = 0;
  Vec<Vec<ZZ_p>> result;
  result.SetLength(type.length());
  for (long i = 0; i < type.length(); i++){
    Vec<ZZ_p> f;
    f.SetLength(type[i]+1);
    for (long j = 0; j < type[i]+1; j++)
      f[j] = v[acc+j];
    acc += type[i] + 1;
    result[i] = f;
  }
  return result;
}

/*----------------------------------------------------------------*/
/* splits v according to the current type                         */
/*----------------------------------------------------------------*/
Vec<zz_pX> hermite_pade::split_on_type(const Vec<zz_p> &v){
  long acc = 0;
  Vec<zz_pX> result;
  for (long i = 0; i < type.length(); i++){
    Vec<zz_p> f;
    for (long j = 0; j < type[i]+1; j++)
      f.append(v[acc+j]);
    acc += type[i] + 1;
    result.append(conv<zz_pX>(f));
  }
  return result;
}

/*----------------------------------------------------------------*/
/* splits v according to the current type                         */
/*----------------------------------------------------------------*/
Vec<ZZX> hermite_pade::split_on_type(const Vec<ZZ> &v){
  long acc = 0;
  Vec<ZZX> result;
  for (long i = 0; i < type.length(); i++){
    Vec<ZZ> f;
    for (long j = 0; j < type[i]+1; j++)
      f.append(v[acc+j]);
    acc += type[i] + 1;
    result.append(conv<ZZX>(f));
  }
  return result;
}

/*----------------------------------------------------------------*/
/* switches ZZ_p to ZZ mod p^(2^n)                                */
/* does not save the current context                              */
/*----------------------------------------------------------------*/
void hermite_pade::switch_context(long n){
  if (n < vec_M.length()){ // it has already been computed
    ZZ_p::init(p_powers[n]);
    level = n;
  } 
  else{
    if (n != vec_M.length())
      throw "length mismatch for table of bmc";

    // calculating the new power of p
    ZZ p_new(p);
    long pow2 = 1L << n;
    power(p_new, p_new, pow2); // 2^n isn't going to be very large
    p_powers.append(p_new);
    ZZ_p::init(p_new);

    // cout << "before making bmc" << endl;
    // string s;
    // getline(cin,s);
    // creating the new bivariate modular comp
    set_up_bmc();
    // cout << "after making bmc" << endl;

    // computing w mod p^(2^n)
    ZZ new_w;
    lift_root_of_unity(new_w, this->w, order, p, pow2);
    vec_w.append(new_w);
    ZZ_p w_p, c_p, d_p;
    conv(w_p, new_w);
    conv(c_p, c);
    conv(d_p, d);


    ZZ_pX_Multipoint_FFT X_new (w_p, c_p, sizeX);
    ZZ_pX_Multipoint_FFT Y_new (w_p, d_p, sizeY);
    // update
    
    vec_X_int.append(X_new);
    vec_Y_int.append(Y_new);
    level = n;

  }
}


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
void hermite_pade::generators_cauchy(Mat<ZZ_p> & X, Mat<ZZ_p> & Y){

  Mat<ZZ_p> G_ZZ_p, H_ZZ_p;
  Vec<ZZ_p> tmp_v, tmp_w;
  generators_hankel_last_row_last_column(G_ZZ_p, H_ZZ_p, tmp_w, tmp_v);

  long n = G_ZZ_p.NumRows();
  long m = H_ZZ_p.NumRows();
  long alpha = G_ZZ_p.NumCols();

  Vec<ZZ_p> f_ZZ_p = conv<Vec<ZZ_p>>(this->f);
  Vec<ZZ_p> e_ZZ_p = conv<Vec<ZZ_p>>(this->e);

  ZZ_p w_ZZ_p, d_ZZ_p;
  w_ZZ_p = conv<ZZ_p>(vec_w[level]);
  d_ZZ_p = conv<ZZ_p>(d);

  X.SetDims(NumRows(), alpha+2);
  Y.SetDims(NumCols(), alpha+2);


  Vec<ZZ_p> tmp_v2;
  for (long j = 0; j < alpha; j++){
    ZZ_pX tmp_p;
    tmp_p.rep.SetLength(n);
    ZZ_p* coef_p = tmp_p.rep.elts();
    for (long i = 0; i < n; i++)
      coef_p[i] = G_ZZ_p[i][j];
    tmp_p.normalize();
    vec_X_int[level].evaluate(tmp_v2, tmp_p);
    for (long i = 0; i < n; i++)
      X[i][j] = tmp_v2[i] * e_ZZ_p[i];
  }

  ZZ_p tmp_z = to_ZZ_p(1);
  for (long i = 0; i < n; i++){
    X[i][alpha] = e_ZZ_p[i] * (power(tmp_z, n)-1);
    tmp_z = tmp_z * w_ZZ_p;
  }


  ZZ_pX tmp_p;
  tmp_p.rep.SetLength(n);
  ZZ_p* coef_p = tmp_p.rep.elts();
  for (long i = 0; i < n; i++)
    coef_p[i] = tmp_v[i];
  tmp_p.normalize();
  vec_X_int[level].evaluate(tmp_v2, tmp_p);
  for (long i = 0; i < n; i++)
    X[i][alpha+1] = tmp_v2[i] * e_ZZ_p[i];

  Vec<ZZ_p> tmp_w2;
  for (long j = 0; j < alpha; j++){
    ZZ_pX tmp_q;
    tmp_q.rep.SetLength(m);
    ZZ_p* coef_q = tmp_q.rep.elts();
    for (long i = 0; i < m; i++)
      coef_q[i] = H_ZZ_p[i][j];
    tmp_q.normalize();
    vec_Y_int[level].evaluate(tmp_w2, tmp_q);
    for (long i = 0; i < m; i++)
      Y[i][j] = tmp_w2[i] * f_ZZ_p[i];
  }
  
  ZZ_pX tmp_q;
  tmp_q.rep.SetLength(m);
  ZZ_p* coef_q = tmp_q.rep.elts();
  for (long i = 0; i < m; i++)
    coef_q[i] = tmp_w[i];
  tmp_q.normalize();
  vec_Y_int[level].evaluate(tmp_w, tmp_q);
  for (long i = 0; i < m; i++)
    Y[i][alpha] = tmp_w[i] * f_ZZ_p[i];

  tmp_z = d_ZZ_p;
  for (long i = 0; i < m; i++){
    Y[i][alpha+1] = -f_ZZ_p[i]*power(tmp_z, m);
    tmp_z = tmp_z * w_ZZ_p;
  }

}


/*----------------------------------------------------------------*/
/* calls create_bmc and appends the result to vec_M               */
/*----------------------------------------------------------------*/
void hermite_pade::set_up_bmc(){
  vec_M.append(create_bmc());
}

/*----------------------------------------------------------------*/
/* multiplies b by the matrix Y_int^t D_f                         */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_Y_right(const Vec<ZZ_p> &b){
  Vec<ZZ_p> x;
  Vec<ZZ_p> f = conv<Vec<ZZ_p>>(this->f); // TODO: e and f should be called e_ZZ and f_ZZ?
  mul_diagonal(x, f, b);
  vec_Y_int[level].mul_left(x, x);
  return x;
}
/*----------------------------------------------------------------*/
/* multiplies b by the matrix D_f Y_int                           */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_Y_left(const Vec<ZZ_p> &b){
  Vec<ZZ_p> b1;
  vec_Y_int[level].mul_right(b1, b);
  Vec<ZZ_p> f = conv<Vec<ZZ_p>>(this->f);
  Vec<ZZ_p> x;
  mul_diagonal(x, f, b1);
  return x;
}


/*----------------------------------------------------------------*/
/* multiplies b by the matrix D_e X_int                           */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_X_right(const Vec<ZZ_p> &b){
  Vec<ZZ_p> b1;
  vec_X_int[level].mul_right(b1, b);
  Vec<ZZ_p> e = conv<Vec<ZZ_p>>(this->e);
  Vec<ZZ_p> x;
  mul_diagonal(x, e, b1);
  return x;
}
/*----------------------------------------------------------------*/
/* multiplies b by the matrix X_int^t D_e                         */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_X_left(const Vec<ZZ_p> &b){
  Vec<ZZ_p> x;
  Vec<ZZ_p> e = conv<Vec<ZZ_p>>(this->e); // TODO: e and f should be called e_ZZ and f_ZZ?
  mul_diagonal(x, e, b);
  vec_X_int[level].mul_left(x, x);
  return x;
}

/*----------------------------------------------------------------*/
/* multiplies b by the matrix M                                   */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_M_right(const Vec<ZZ_p> &b){
  Vec<ZZ_p> flipped = flip_on_type(b);
  return vec_M[level]->mul_right(flipped);
}
/*----------------------------------------------------------------*/
/* left-multiplies b by the matrix M                              */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_M_left(const Vec<ZZ_p> &b){
  Vec<ZZ_p> flipped = vec_M[level]->mul_left(b);
  return flip_on_type(flipped);
}


/*----------------------------------------------------------------*/
/* multiplies b by the matrix CL =  D_e X_int M Y_int^t D_f       */
/* (CL is cauchy-geometric-like)                                  */
/* b need not have size CL.NumCols()                              */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mulA_right(const Vec<ZZ_p>& b){
  double t = GetTime();
  double t2;
  Vec<ZZ_p> b_loc = b;
  b_loc.SetLength(sizeY, ZZ_p(0)); // pad it
  t2 = GetTime();
  Vec<ZZ_p> x = mul_Y_right(b_loc);
  time_mul_y_right += GetTime() - t2;
  t2 = GetTime();
  x = mul_M_right(x);
  time_mul_m_right += GetTime() - t2;
  t2 = GetTime();
  x = mul_X_right(x);
  time_mul_x_right += GetTime() - t2;
  time_mulA += GetTime() - t;
  return x;
}
/*----------------------------------------------------------------*/
/* multiplies b by the matrix M                                   */
/*----------------------------------------------------------------*/
Mat<ZZ_p> hermite_pade::mulA_right(const Mat<ZZ_p> &b){
  Mat<ZZ_p> output;
  output.SetDims(NumRows(), b.NumCols());
  Vec<ZZ_p> inv, outv;
  inv.SetLength(b.NumRows());
  for(long i = 0; i < b.NumCols(); i++){
    for (long j = 0; j < b.NumRows(); j++)
      inv[j] = b[j][i];
    outv = mulA_right(inv);
    for (long j = 0; j < NumRows(); j++)
      output[j][i] = outv[j];
  }
  return output;
}

/*----------------------------------------------------------------*/
/* left-multiplies b by the matrix CL =  D_e X_int M Y_int^t D_f  */
/* (CL is cauchy-geometric-like)                                  */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mulA_left(const Vec<ZZ_p>& b){
  Vec<ZZ_p> b1=b;
  b1.SetLength(e.length());
  double time = GetTime();
  Vec<ZZ_p> x = mul_X_left(b1);
  time_mul_x_left += GetTime() - time;
  time = GetTime();
  x = mul_M_left(x);
  time_mul_m_left += GetTime() - time;
  time = GetTime();
  x = mul_Y_left(x);
  time_mul_y_left += GetTime() - time;

  return x;
}
/*----------------------------------------------------------------*/
/* left-multiplies b by the matrix M                              */
/*----------------------------------------------------------------*/
Mat<ZZ_p> hermite_pade::mulA_left(const Mat<ZZ_p> &b){
  Mat<ZZ_p> output;
  output.SetDims(NumCols(), b.NumCols());
  Vec<ZZ_p> inv, outv;
  inv.SetLength(b.NumRows());

  for(long i = 0; i < b.NumCols(); i++){
    for (long j = 0; j < b.NumRows(); j++)
      inv[j] = b[j][i];
    outv = mulA_left(inv);
    for (long j = 0; j < NumCols(); j++)
      output[j][i] = outv[j];
  }
  return output;
}

/*----------------------------------------------------------------*/
/* if Mx = b mod p^(2^{n-1}), updates x so that Mx = b mod p^(2^n)*/
/*----------------------------------------------------------------*/
void hermite_pade::update_solution(Vec<ZZ>& x, const Vec<ZZ_p> &b, long n){

  Vec<ZZ_p> r = mulA_right(conv<Vec<ZZ_p>>(x)) - b; // mod p^(2^n)
  Vec<ZZ> r_ZZ = conv<Vec<ZZ>>(r);
  for (long i = 0; i < r.length(); i++)
    r_ZZ[i] = r_ZZ[i] / p_powers[n-1];

  Vec<ZZ> x_1;
  switch(mode){
  case 0:
    DAC(x_1, r_ZZ, n-1);
    break;
  case 1:
    Dixon(x_1, r_ZZ, n-1);
    break;
  case 2:
    solve_newton(x_1, r_ZZ, n-1);
    break;
  default:
    throw "solving mode not implemented\n";
  }

  x = x - p_powers[n-1] * x_1;
}

/*----------------------------------------------------------------*/
/* solves Mx = b mod p^(2^n) by Newton iteration                  */
/*----------------------------------------------------------------*/
void hermite_pade::solve_newton(Vec<ZZ>& x, const Vec<ZZ> &b, long n){
  while (lift_invA.length() <= n)
    Newton(lift_invA.length());

  long old_n = level;
  switch_context(n);

  Vec<ZZ_p> x_ZZ_p;

  lift_invA[n].mul_right(x_ZZ_p, conv<Vec<ZZ_p>>(b));
  x = conv<Vec<ZZ>>(x_ZZ_p);
  switch_context(old_n);
}

/*----------------------------------------------------------------*/
/* computes the inverse of A (cauchy matrix) mod p^(2^n)          */
/* assumes the inverse mod p^(2^(n-1)) is known                   */
/*----------------------------------------------------------------*/
void hermite_pade::Newton(long n){

  if (n == 0){
    zz_pContext push;
    ctx.restore();
    long old_n = level;
    switch_context(0);
    lift_invA.SetLength(0);
    ZZ_p_cauchy_like_geometric C(conv<Mat<ZZ_p>>(conv<Mat<ZZ>>(invA.G)), conv<Mat<ZZ_p>>(conv<Mat<ZZ>>(invA.H)), 
				 conv<ZZ_p>(conv<ZZ>(invA.C.u1)), conv<ZZ_p>(conv<ZZ>(invA.C.v1)), conv<ZZ_p>(conv<ZZ>(invA.C.rho)) );
    lift_invA.append(C);
    switch_context(old_n);
    return;
  }
  // set up the first inverse mod p
  long old_n = level;

  switch_context(n);

  // gens of cauchy matrix taken mod p^(2^n)
  Mat<ZZ_p> X, Y, X1, Y1, X2, Y2, X3, Y3;
  generators_cauchy(X, Y);


  ZZ_p_cauchy_like_geometric C_n(conv<Mat<ZZ_p>>(conv<Mat<ZZ>>(lift_invA[n-1].G)), 
  				 conv<Mat<ZZ_p>>(conv<Mat<ZZ>>(lift_invA[n-1].H)), 
				 conv<ZZ_p>(d), conv<ZZ_p>(c), conv<ZZ_p>(vec_w[n]) );

  // cout << "cauchy dimensions " << C_n.NumRows() << " " << C_n.NumCols() << " " << C_n.NumGens() << endl;
  // cout << "X dimensions " << X.NumRows() << " " << X.NumCols() << " " << endl;
  double time = GetTime();
  C_n.mul_right(X1, X);
  time_cauchy_right += GetTime() - time;
  X2 = mulA_right(X1);
  time = GetTime();
  C_n.mul_right(X3, X2);
  time_cauchy_right += GetTime() - time;
  time = GetTime();
  C_n.mul_left(Y1, Y);
  time_cauchy_left += GetTime() - time;
  Y2 = mulA_left(Y1);
  time = GetTime();
  C_n.mul_left(Y3, Y2);
  time_cauchy_left += GetTime() - time;

  lift_invA.append(  ZZ_p_cauchy_like_geometric(-(2*X1-X3), (2*Y1-Y3), conv<ZZ_p>(d), conv<ZZ_p>(c), conv<ZZ_p>(vec_w[n])) );

#if false 
  Mat<ZZ> Hd;
  to_dense(Hd, MH_ZZ);
  Mat<ZZ_p> Xi, Yi, De, Df, H, C, invC;
  vec_X_int[n].to_dense(Xi);
  vec_Y_int[n].to_dense(Yi);
  De.SetDims(e.length(), e.length());
  for (long i = 0; i < e.length(); i++)
    De[i][i] = conv<ZZ_p>(e[i]);
  Df.SetDims(f.length(), f.length());
  for (long i = 0; i < f.length(); i++)
    Df[i][i] = conv<ZZ_p>(f[i]);
  H = conv<Mat<ZZ_p>>(Hd);
  C = De*Xi * H * transpose(Df*Yi);
  lift_invA[n].to_dense(invC);
  cout << invC*C << endl;
#endif

  switch_context(old_n);
}

/*----------------------------------------------------------------*/
/* solves for Mx = b mod p^(2^n)                                  */
/*----------------------------------------------------------------*/
void hermite_pade::DAC(Vec<ZZ> &x, const Vec<ZZ>& b_in, long n){
  if (n == 0){ // we are mod p
    zz_pContext push;
    ctx.restore();
    Vec<zz_p> x_zz_p;
    Vec<zz_p> b_zz_p = conv<Vec<zz_p>>(b_in);
    b_zz_p.SetLength(rank); // TODO: why do we need this?
    invA.mul_right(x_zz_p, b_zz_p);
    x = conv<Vec<ZZ>>(x_zz_p);
    return;
  }

  long old_n = level;
  switch_context(n); 

  Vec<ZZ_p> b_n = conv<Vec<ZZ_p>>(b_in);
  Vec<ZZ> b = conv<Vec<ZZ>>(b_n); // this reduces b_in mod p^(2^n).
  DAC(x, b, n-1); 

  update_solution(x, b_n, n);
  switch_context(old_n); 
}

/*----------------------------------------------------------------*/
/* solves for Mx = b mod p^(2^n) using Dixon's algorithm          */
/*----------------------------------------------------------------*/
void hermite_pade::Dixon (Vec<ZZ> &x, Vec<ZZ> &b_in, long n){
  // Vector of the x's
  Vec<Vec<ZZ>> vec_x;
  
  // computing x_0
  zz_pContext push;
  ctx.restore();
  Vec<zz_p> x_zz_p;
  Vec<zz_p> b_zz_p = conv<Vec<zz_p>>(b_in);
  
  invA.mul_right(x_zz_p, b_zz_p);
  vec_x.append(conv<Vec<ZZ>>(x_zz_p));
  
  auto old_n = level;
  long t = power_long(2, n);
  switch_context(n);
  
  for (long i = 1; i < t; ++i){
    Vec<ZZ_p> b_ZZ_p = conv<Vec<ZZ_p>>(b_in);
    auto c_ZZ_p = mulA_right(conv<Vec<ZZ_p>>(vec_x[i-1]));
    b_ZZ_p = b_ZZ_p - c_ZZ_p;
    b_in = conv<Vec<ZZ>>(b_ZZ_p);
    for (long j = 0; j < b_ZZ_p.length(); j++)
      b_in[j] = b_in[j] / p;
    b_zz_p = conv<Vec<zz_p>>(b_in);
    invA.mul_right(x_zz_p, b_zz_p);
    vec_x.append(conv<Vec<ZZ>>(x_zz_p));
  }
  
  switch_context(n);
  Vec <ZZ_p>x_ZZ_p = conv<Vec<ZZ_p>>(vec_x[0]); // running total
  ZZ p_running = ZZ(p);
  for (long i = 1; i < t; i++){
    Vec<ZZ_p> temp = conv<Vec<ZZ_p>>(vec_x[i]);
    x_ZZ_p += conv<ZZ_p>(p_running) * temp;
    p_running *= p;
  }
  x = conv<Vec<ZZ>>(x_ZZ_p);
  switch_context(old_n);
}

/*----------------------------------------------------------------*/
/* checks if every entry can be reconstructed                     */
/* if so, check whether the solution cancels the system mod p2    */
/*----------------------------------------------------------------*/
bool hermite_pade::reconstruct_and_check(Vec<ZZX> & sol_poly, const Vec<ZZ_p> &v, long n){
  double t = GetTime();
  Vec<Vec<ZZ>> sol;
  sol.SetLength(0);

  Vec<long> flip_wit = flip_on_type(witness); // local variable?
  
  ZZ bound;
  if (n == 0)
    bound = SqrRoot(ZZ(p));
  else
    bound = p_powers[n-1];

  ZZ_p denom = to_ZZ_p(0);
  for (long i = 0; i < v.length(); i++)
    if (v[i] != 0){
      denom = v[i];
      break;
    }
  ZZ_p i_den = 1 / denom;

  for (long i = 0; i < v.length(); i++){
    ZZ a,b;
    try{
      long result = ReconstructRational(a, b, conv<ZZ>(v[i] * i_den), p_powers[n], bound, bound);
      if (result == 0){
        time_recon += GetTime() - t;
        time_recon_all += GetTime()-t; 
	return false;
      }
      Vec<ZZ> temp;
      temp.append(a);
      temp.append(b);
      sol.append(temp);

      if (stop_criterion == 2) {
	ZZ_pPush push;
	ctx2.restore();
	if (conv<long>(conv<ZZ>(conv<ZZ_p>(a) / conv<ZZ_p>(b))) != flip_wit[i]){
	  cout << "early abort in ratrecon with n=" << n << " i=" << i << endl;
	  time_recon_all += GetTime()-t;
	  return false;
	}
      }

    }
    catch(...){
      time_recon += GetTime() - t;
      time_recon_all += GetTime()-t;
      return false;
    }
  }
  
  time_recon += GetTime() - t;
  sol = flip_on_type(sol);

  // 0 = re-evaluate mod another prime
  // 1 = twice the same result
  // 2 = compare to witness

  // for stop_criterion = 1
  static Vec<Vec<ZZ>> sol_old;

  switch (stop_criterion){
  case 0:
    {
      double t2 = GetTime();
      ZZ_pContext push;
      ctx2.restore();
      
      Vec<ZZ_p> sol_ZZ_p;
      for (long int i = 0; i < sol.length(); i++)
	sol_ZZ_p.append(conv<ZZ_p>(sol[i][0]) / conv<ZZ_p>(sol[i][1]));
      
      auto bmc = create_bmc();
      auto x = bmc->mul_right(sol_ZZ_p);
      time_check_p2 += GetTime() - t2;
      for (long int i = 0; i < x.length(); i++)
	if (x[i] != ZZ_p(0)){
	  cout << "failed ratrecon with n=" << n << endl;
	  time_recon_all += GetTime()-t;
	  return false;
	}
    }
    break;
  case 1:
    {
      if (sol != sol_old){
	sol_old = sol;
	cout << "failed ratrecon with n=" << n << endl;
	time_recon_all += GetTime()-t;
	return false;
      }
    }
    break;
  case 2:
    {
      ZZ_pPush push;
      ctx2.restore();
      Vec<ZZ_p> sol_ZZ_p;
      for (long int i = 0; i < sol.length(); i++)
 	sol_ZZ_p.append(conv<ZZ_p>(sol[i][0]) / conv<ZZ_p>(sol[i][1]));
      if ( conv<Vec<long>>(conv<Vec<ZZ>>(sol_ZZ_p)) != witness ){
	cout << "failed ratrecon using witness with n=" << n << endl;
	time_recon_all += GetTime()-t;
	return false;
      }
    }
    break;
  }

  double time = GetTime();
  Vec<ZZ> sol_ZZ;
  sol_ZZ.SetLength(sol.length());
  ZZ ell{1};
  for (long i = 0; i < sol.length(); i++)
    ell = (sol[i][1] * ell) / GCD(sol[i][1], ell);
  for (long i = 0; i < sol.length(); i++)
    sol_ZZ[i] = (sol[i][0] * ell) / sol[i][1]; 
  div_in_recon += GetTime()-time;   

  sol_poly = split_on_type(sol_ZZ);
  time_recon_all += GetTime()-t;
  return true;
}

/*----------------------------------------------------------------*/
/* computes a random solution of the system modulo p              */
/*----------------------------------------------------------------*/
void hermite_pade::random_solution_mod_p(Vec<zz_pX> & pol){
  pol.SetLength(type.length());

  Vec<ZZ_p> ex1; // mult with A to get a column
  Vec<ZZ_p> add1;
  Vec<ZZ> sol1;

  ex1.SetLength(sizeY);
  add1.SetLength(sizeY - rank);
  for (long i = 0; i < rank; i++)
    ex1[i] = 0;
  for (long i = 0; i < add1.length(); i++){
    add1[i] = random_ZZ_p();
    ex1[rank+i] = add1[i];
  }
  ex1 = mulA_right(ex1);

  Vec<ZZ> tmp_in = conv<Vec<ZZ>>(ex1);

  DAC(sol1, tmp_in, 0);

  conv(ex1, sol1);
  ex1.SetLength(sizeY, ZZ_p(0));
  
  for (long i = 0; i < add1.length(); i++)
    ex1[rank+i] = -add1[i];
	
  ex1 = flip_on_type(mul_Y_right(ex1));
  Vec<zz_p> v = conv<Vec<zz_p>>(conv<Vec<ZZ>>(ex1));

  zz_p denom = to_zz_p(0);
  for (long i = 0; i < v.length(); i++)
    if (v[i] != 0){
      denom = v[i];
      break;
    }
  if (denom == 0)
    return;
  for (long i = 0; i < v.length(); i++)
    v[i] /= denom;
  pol = split_on_type(v);
}

/*----------------------------------------------------------------*/
/* computes a solution of the system modulo p                     */
/* assumes that the nullity is one                                */
/*----------------------------------------------------------------*/
void hermite_pade::random_solution(Vec<ZZX> &sol_poly){
	double time_all = GetTime();
  if ((NumCols() - Rank()) != 1)
    throw "lifting works only in case of nullity 1";

  zz_pPush zz_p_push;
  ZZ_pPush ZZ_p_push;

  ctx.restore();
  long n = 0; // start at p^2^n
  switch_context(n);


  Vec<ZZ_p> extractor; // mult with A to get a column
  extractor.SetLength(sizeY, ZZ_p(0));
  extractor[rank] = 1; // just for now, take the last column

  Vec<ZZ_p> b = mulA_right(extractor); // b is the last column of A
  Vec<ZZ> x_ZZ, b_ZZ;
  b_ZZ = conv<Vec<ZZ>> (b);

  DAC(x_ZZ, b_ZZ, n); // solution mod p

  Vec<ZZ_p> x = conv<Vec<ZZ_p>> (x_ZZ);
  x.SetLength(sizeY, ZZ_p(0));
  x[rank] = -1;
  Vec<ZZ_p> soln = mul_Y_right(x);
	
	double time_loop = GetTime();
  // loop until we get enough prec
  while (!reconstruct_and_check(sol_poly, soln, n)){
    switch_context(++n);

    Vec<ZZ_p> extractor; // mult with A to get a column
    extractor.SetLength(sizeY, ZZ_p(0));
    extractor[rank] = 1; // just for now, take the last column
    Vec<ZZ_p> b = mulA_right(extractor); // b is the last column of A

    update_solution(x_ZZ, b, n);

    Vec<ZZ_p> x = conv<Vec<ZZ_p>> (x_ZZ);
    x.SetLength(sizeY,ZZ_p(0));
    x[rank] = -1;
    soln.kill();
    soln = mul_Y_right(x);

    ZZ_p first;
    for (long i = 0; i < soln.length(); i++)
      if (soln[i]._ZZ_p__rep % p_powers[0] != ZZ(0)){
        first = 1/soln[i];
        break;
      }
    soln *= first;
  }
  cout << "everything: " << GetTime() - time_all << endl;
  cout << "loop: " << GetTime() - time_loop << endl;
  cout << "mul time: " << time_mulA << endl;
  if (time_spent_lower != 0)
    cout << "Time for lower: " << time_spent_lower << endl;
  cout << "mul M right: " << time_mul_m_right << endl;
  cout << "mul X right: " << time_mul_x_right << endl;
  cout << "mul Y right: " << time_mul_y_right << endl;
  if (mode == 2){
    cout << "mul M left: " << time_mul_m_left << endl;
    cout << "mul X left: " << time_mul_x_left << endl;
    cout << "mul Y left: " << time_mul_y_left << endl;
    cout << "time cauchy left: " << time_cauchy_left << endl;
    cout << "time cauchy right: " << time_cauchy_right << endl;
  }
  cout << "total reconstruction time: " << time_recon_all << endl;
  cout << "reconstruction only: " << time_recon << endl;
  cout << "divisions in recon: " << div_in_recon << endl; 
  cout << "checking p2: " << time_check_p2 << endl;
}


/*----------------------------------------------------------------*/
/* Rank                                                           */
/*----------------------------------------------------------------*/
long hermite_pade::Rank() const {
  return rank;
}

/*----------------------------------------------------------------*/
/* Dimensions: number of rows                                     */
/*----------------------------------------------------------------*/
long hermite_pade::NumRows() const {
  return sizeX;
}

/*----------------------------------------------------------------*/
/* Dimensions: number of columns                                  */
/*----------------------------------------------------------------*/
long hermite_pade::NumCols() const {
  return sizeY;
}

/*----------------------------------------------------------------*/
/* sets up the field, contexts, ..                                */
/*----------------------------------------------------------------*/
hermite_pade::hermite_pade(long fft_index) {
  level = 0;
  zz_p::FFTInit(fft_index+1);
  ctx2 = ZZ_pContext(ZZ(zz_p::modulus())); 
  zz_p::FFTInit(fft_index);
  ctx = zz_pContext(INIT_FFT, fft_index);
  p = zz_p::modulus();
  ZZ_p::init(ZZ(p));
  p_powers.append(ZZ(p));
}


/*----------------------------------------------------------------*/
/* mode switch. 0-DAC, 1-Dixon, 2-Newton                          */
/*----------------------------------------------------------------*/
void hermite_pade::switch_mode(long i){
  if (i < 0 || i > 2) 
    throw "mode not supported";
  mode = i;
}


/*----------------------------------------------------------------*/
/* does nothing                                                   */
/*----------------------------------------------------------------*/
hermite_pade::hermite_pade(){
}

#ifndef LIN_SEQ
#define LIN_SEQ

#include <NTL/lzz_pX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/vec_lzz_p.h>
#include <cmath>
#include "lzz_pX_CRT.h"
#include <NTL/ZZ.h>
#include "lzz_pXY.h"
#include <map>

NTL_CLIENT

// class to compute a sequence (a_i)
class lin_seq{
	long p; // prime field
  zz_pX P; // annihiliating polynomial
	vec_pair_zz_pX_long sq_free_pols; // Q_i's	
	zz_pX_CRT* crt;
	Vec<zz_p> init; // vector of initial conditions
	Vec<zz_pX> expanded; // stores Qi^mi
	Vec<zz_pXModulus> sq_free_mods; // stores zz_pXModulus of Q_i's

	// precompute some quantities
	Vec<zz_p> factorials;
	Vec<zz_p> inv_factorials;
	
/** DIVISION ROUTINES ******************************/
	// compute x^D mod Q_i^m_i
	void mod(zz_pX&, const ZZ &D, const long i);
	
	// generates m binomial coeffs (D choose i)
	void gen_binomial(vec_zz_p&, const long m, const ZZ &D);
	
	// repeatedly squares p to k
	void repeated_sq(zz_pX &, const zz_pX& p, const ZZ &k);
	
	// repeatedly squares p mod h to k
	void repeated_sq_mod(zz_pX &, const zz_pX& p, const zz_pX& h, const ZZ &k);
	
	// returns the binary representation of n
	Vec<int> get_binary(const ZZ &n); 
	
	// compues x^D mod (x-1)^m
	void compute_x(zz_pX &, const ZZ &D, const long m);

/*****Shifting **********************************************/
	// given P, computes P(x+1) using the matrix method
	void shift(zz_pX &result, const zz_pX &P);

	// given P and s, computes P(x+s)
	void shift_by(zz_pX &result, const zz_pX &P, const zz_p &s);

	// computes the shift of v
	void mult_shift(Vec<zz_p> &result, const Vec<zz_p> &v);
	
	// computes the shift of f = a_0 y^t+a_1 y^{t-1} x^2+...+a_m y^{t-m} x^m
	// (i.e. f is ordered by powers of x)
	void shift_by_alpha (zz_pXY &result, const zz_pXY &f, const long ind, const int dir);
/************************************************************/

/*****Transpose Shifting*************************************/
	// given P, computes P(x+1) using the matrix method
	void T_shift(zz_pX &result, const zz_pX &P);

	// given P and s, computes P(x+s)
	void T_shift_by(zz_pX &result, const zz_pX &P, const zz_p &s);

	// computes the shift of v
	void T_mult_shift(Vec<zz_p> &result, const Vec<zz_p> &v);

  void T_shift_by_alpha (zz_pXY &result, const zz_pXY &f, 
	                       const long ind, const int dir);	
/************************************************************/

	// given the dimension of the shift matrix, computes
	// the polynomial representation of the Hankel matrix
	void get_h_p (zz_pX &result, const long &n);
	
	// computes taylor expansion of a at p, where deg(p) = 1
	void t_expand(Vec<zz_p>&, const zz_pX& a, const zz_pX& p, const ZZ &k);
	
	// computes x^D mod (h(y), (x-y)^m)
	void bivariate_mod (Vec<zz_pX> &, const ZZ &D, const zz_pX& h, const long m);
	
	// computes the square free decomposition
	void sqfree_decomp(vec_pair_zz_pX_long& result,const zz_pX f);

	void untangling_transpose(zz_pX &result, const zz_pXY &in, const long i); 

	// untangling routines
	zz_pX untangling_transpose(zz_pXY u, const map<long,zz_pX>& hbar_powers, long m); 
/*****************************************************/

public:
	lin_seq(const zz_pX& P, const Vec<zz_p>& init, const long &p);
	
	// returns the D-th entry specified by P
	void get_entry(zz_p&, const ZZ &D);
	// naively computes D
	void get_entry_naive(zz_p&, const ZZ &D);

	~lin_seq(){delete crt;}

	// solves for a coefficients (a_0,...,a_t) st
	// a_0*N + a_1*x*N + ... + a_t*x^t*N = M mod P
	// where N/P generates lhs, and M/P generates rhs
	// NOTE: P cancels lhs and rhs, so N/P and M/P is
	// a series over 1/T
	void hankel_solve(zz_pX &result, const Vec<zz_p> &lhs, 
	                  const Vec<zz_p> &rhs, const long ind);  

};


#endif

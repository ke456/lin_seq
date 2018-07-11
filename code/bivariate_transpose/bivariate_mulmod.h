#ifndef __BIVARIATE_MULMOD_H
#define __BIVARIATE_MULMOD_H

#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <string.h>

NTL_CLIENT

// initializes polynomials in sage
void init_sage();

// prints a univariate polynomial in y
void print(const zz_pX& A);

// prints a polynomial in x, with coefficients in y
void print(const Vec<zz_pX>& A);

// prints the assignment of variable name to polynomial A (in sage)
//void sage_assign(const zz_pX& A, const string& name);

// prints the assignment of variable name to bivariate polynomial A (in sage)
void sage_assign(const Vec<zz_pX>& A, const string& name);

// random vector of length m with coefficients of degree < d
void random(Vec<zz_pX>& res, long d, long m);

// multiplies by the matrix Yp = [h1 h2 ... hm]
//                               [h2 ... hm 0 ]
//                               [...         ]
//                               [hm          ]
void mul_Yp(zz_pX& res, const zz_pX& a, const zz_pX& h);

// multiplies by the inverse of Yp = [h1 h2 ... hm]
//                                   [h2 ... hm 0 ]
//                                   [...         ]
//                                   [hm          ]
// assumes we are given 1/rev(h) mod x^m
void mul_inv_Yp(zz_pX& res, const zz_pX& a, const zz_pX& h, const zz_pX& inv_rev_h);

// a and b have length EXACTLY m
// multiplies a by b modulo h(y), x^m
// result has length m
void bivariate_product(Vec<zz_pX>& res, const Vec<zz_pX>& a, const Vec<zz_pX>& b, const zz_pX& h);

// same, using kronecker substitution
void bivariate_product_kronecker(Vec<zz_pX>& c, const Vec<zz_pX>& a, const Vec<zz_pX>& b, const zz_pX& h);


// a and b have length EXACTLY m
// res = a^t b
// result has length m
void bivariate_transposed_product(Vec<zz_pX>& res, const Vec<zz_pX>& a, const Vec<zz_pX>& b, const zz_pX& h);

#endif

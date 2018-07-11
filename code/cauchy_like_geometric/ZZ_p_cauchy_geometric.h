#ifndef __ZZ_P_CAUCHY_GEOMETRIC_H
#define __ZZ_P_CAUCHY_GEOMETRIC_H

#include <NTL/ZZ_pX.h>

#include "ZZ_p_toeplitz.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Cauchy matrices on geometric progressions          */
/*----------------------------------------------------*/
/*----------------------------------------------------*/
class ZZ_p_cauchy_geometric{
 private:
  long m, n;
 
 public:
  ZZ_p u1, v1, rho, sqrt_rho;
  Vec<ZZ_p> vec_toeplitz;
  ZZ_p_toeplitz t;
  Vec<ZZ_p> powers_irho;

  /*---------------------------------------------------*/
  /* default constructor                               */
  /*---------------------------------------------------*/
  ZZ_p_cauchy_geometric();

  /*---------------------------------------------------*/
  /* constructor                                       */
  /*---------------------------------------------------*/
  ZZ_p_cauchy_geometric(const ZZ_p& a1, const ZZ_p& b1, const ZZ_p& rho, long mm, long nn);

  /*---------------------------------------------------*/
  /* dimensions                                        */
  /*---------------------------------------------------*/
  long NumRows() const;
  long NumCols() const;

  /*---------------------------------------------------*/
  /* computes output = M*input                         */
  /*---------------------------------------------------*/
  void mul_right(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const;
  void mul_right(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const;

  /*---------------------------------------------------*/
  /* computes output = M*input, without the diagonal   */
  /*---------------------------------------------------*/
  void mul_right_simple(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const;
  void mul_right_simple(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const;

  /*---------------------------------------------------*/
  /* computes output = M^t*input                       */
  /*---------------------------------------------------*/
  void mul_left(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const;
  void mul_left(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const;

  /*---------------------------------------------------*/
  /* computes output = M^t*input, without the diagonal */
  /*---------------------------------------------------*/
  void mul_left_simple(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const;
  void mul_left_simple(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const;

  /*---------------------------------------------------*/
  /* M as a dense matrix                               */
  /*---------------------------------------------------*/
  void to_dense(Mat<ZZ_p>& M) const;

};


/*---------------------------------------------------*/
/* computes                                          */
/* 1/(u1-v1 rho^(-m+1)) ... 1/(u1-v1 rho^(n-1))      */
/* these are the entries of the toeplitz matrix      */
/* (with m rows and n columns)                       */
/*---------------------------------------------------*/
void prepare_inverses_cauchy(Vec<ZZ_p>& inverses, const ZZ_p& u1, const ZZ_p& v1, const ZZ_p& rho, long m, long n);



/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Cauchy-like matrices on geometric progressions     */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

class ZZ_p_cauchy_like_geometric{
public:
  ZZ_p_cauchy_geometric C;
  Mat<ZZ_p> G, H;

  /*---------------------------------------------------*/
  /* default constructor                               */
  /*---------------------------------------------------*/
  ZZ_p_cauchy_like_geometric();

  /*---------------------------------------------------*/
  /* constructor                                       */
  /*---------------------------------------------------*/
  ZZ_p_cauchy_like_geometric(const Mat<ZZ_p>& U, const Mat<ZZ_p>& V, const ZZ_p& a1, const ZZ_p& b1, const ZZ_p& rho);

  /*---------------------------------------------------*/
  /* dimensions                                        */
  /*---------------------------------------------------*/
  long NumRows() const;
  long NumCols() const;
  long NumGens() const;

  /*---------------------------------------------------*/
  /* computes output = M*input                         */
  /*---------------------------------------------------*/
  void mul_right(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const;
  void mul_right(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const;

  /*---------------------------------------------------*/
  /* computes output = M^t*input                       */
  /*---------------------------------------------------*/
  void mul_left(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const;
  void mul_left(Mat<ZZ_p>& output, const Mat<ZZ_p>& input) const;

  /*---------------------------------------------------*/
  /* M as a dense matrix                               */
  /*---------------------------------------------------*/
  void to_dense(Mat<ZZ_p>& M) const;
};



#endif

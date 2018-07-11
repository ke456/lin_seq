#ifndef __ZZ_P_TOEPLITZ_H
#define __ZZ_P_TOEPLITZ_H

#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_p.h>


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

class ZZ_p_toeplitz{
 public:
  // n rows, m columns
  long n, m;
  ZZ_pX data, data_rev;

  /*----------------------------------------------------*/
  /* default constructor                                */
  /*----------------------------------------------------*/
  ZZ_p_toeplitz();

  /*----------------------------------------------------*/
  /* input vector is as showed above                    */
  /*----------------------------------------------------*/
  ZZ_p_toeplitz(const Vec<ZZ_p>& input, long rows, long cols);

  /*---------------------------------------------------*/
  /* dimensions                                        */
  /*---------------------------------------------------*/
  long NumRows() const;
  long NumCols() const;

  /*---------------------------------------------------*/
  /* data access                                       */
  /*---------------------------------------------------*/
  const ZZ_p& operator ()(long i, long j) const;

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
  void to_dense(Mat<ZZ_p>& Mdense) const;
};

#endif

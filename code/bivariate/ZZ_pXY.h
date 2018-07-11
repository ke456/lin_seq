#ifndef ZZPXY__H
#define ZZPXY__H

#include <fstream>
#include <iostream>
#include <string.h>

#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>

#include "ZZXY.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a helper class for bivariate polynomials over ZZ_p         */
/* a ZZ_pXY is simply a vector of ZZ_pX                       */
/* with the convention f = sum_i coeffX[i](X) Y^i             */ 
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

class ZZ_pXY{
 public:
  Vec<ZZ_pX> coeffX;

  /*------------------------------------------------------------*/
  /* total degree                                               */
  /*------------------------------------------------------------*/
  long tdeg() const;

  /*------------------------------------------------------------*/
  /* degree in Y                                                */
  /*------------------------------------------------------------*/
  long degY() const;

  /*------------------------------------------------------------*/
  /* degree in X                                                */
  /*------------------------------------------------------------*/
  long degX() const;

  /*------------------------------------------------------------*/
  /* naive evaluation algorithm to compute F(x,g(x)) mod x^t    */ 
  /*------------------------------------------------------------*/
  void eval(ZZ_pX & val, const ZZ_pX & g, long t);

  /*------------------------------------------------------------*/
  /* finds g such that F(x,g(x)) = 0 mod x^t, g(0) = g0         */ 
  /*------------------------------------------------------------*/
  void series_solution(ZZ_pX & g, const ZZ_p & g0, long t);

  /*------------------------------------------------------------*/
  /* resizes the array of coefficients                          */
  /* to remove the trailing entries that are zero, if any       */
  /*------------------------------------------------------------*/
  void normalize();

  ZZ_pXY(){}

  /*------------------------------------------------------------*/
  /* builds from a vector of ZZ_pX                                */
  /*------------------------------------------------------------*/
  ZZ_pXY(const Vec<ZZ_pX> & coeff){
    coeffX = coeff;
    normalize();
  }

  /*------------------------------------------------------------*/
  /* builds from a file                                         */
  /* format: [[a0,a1,..,aN],...,[x0,x1,..,xM]]                  */
  /* where [a0,a1,..,aN] = coeff(f,Y^0) \in Z/pZ[X], etc.       */
  /*------------------------------------------------------------*/
  void read_from_file(const string& filename);

  ~ZZ_pXY(){}
};

/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
istream& operator>>(istream& s, ZZ_pXY& x);
ostream& operator<<(ostream& s, const ZZ_pXY& a);

/*------------------------------------------------------------*/
/* initializes M<y,x>=QQ[y,x] with lex order y > x            */
/*------------------------------------------------------------*/
void magma_init_bi_QQ();

/*------------------------------------------------------------*/
/* prints a poly in M = QQ[y > x]                             */
/*------------------------------------------------------------*/
void magma_output(const ZZ_pXY & a);

/*------------------------------------------------------------*/
/* assigns a poly in M = QQ[y > x]                            */
/*------------------------------------------------------------*/
void magma_assign(const ZZ_pXY & a, const string & name);

/*------------------------------------------------------------*/
/* derivative in Y                                            */
/*------------------------------------------------------------*/
void diffY(ZZ_pXY & dyF, const ZZ_pXY & F);

/*------------------------------------------------------------*/
/* random poly with deg(F,x) <= dx, deg(F,y) <= dy            */
/*------------------------------------------------------------*/
void random(ZZ_pXY & F, long dx, long dy);

/*------------------------------------------------------------*/
/* reduces modulo p                                           */
/*------------------------------------------------------------*/
void conv(ZZ_pXY & f, const ZZXY & fZ);

#endif

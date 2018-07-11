#include <sstream> 
#include <algorithm> 
#include <NTL/tools.h>
#include <NTL/vector.h>

#include "ZZ_pXY.h"
#include "ZZ_pX_extra.h"
#include "magma_output.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a helper class for bivariate polynomials over ZZ_p         */
/* a ZZ_pXY is simply a vector of ZZ_pX                       */
/* with the convention f = sum_i coeffX[i](X) Y^i             */ 
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* total degree                                               */
/*------------------------------------------------------------*/
long ZZ_pXY::tdeg() const {
  long res = -1;
  for (long i = 0; i < coeffX.length(); i++)
    res = max(res, deg(coeffX[i]) + i);
  return res;
}

/*------------------------------------------------------------*/
/* degree in Y                                                */
/*------------------------------------------------------------*/
long ZZ_pXY::degY() const {
  return coeffX.length()-1;
}

/*------------------------------------------------------------*/
/* degree in X                                                */
/*------------------------------------------------------------*/
long ZZ_pXY::degX() const {
  long d = -1;
  for (long i = 0; i < coeffX.length(); i++)
    d = max(d, deg(coeffX[i]));
  return d;
}

/*------------------------------------------------------------*/
/* naive evaluation algorithm to compute F(x,g(x)) mod x^t    */ 
/*------------------------------------------------------------*/
void ZZ_pXY::eval(ZZ_pX & val, const ZZ_pX & g, long t){
  long d = degY();

  if (d == -1){
    clear(val);
    return;
  }

  ZZ_pX gt = trunc(g, t);
  val = coeffX[d];
  for (long i = d-1; i >= 0; i--){
    val = trunc(val * gt, t);
    val = val + trunc(coeffX[i], t);
  }
}

/*------------------------------------------------------------*/
/* finds g such that F(x,g(x)) = 0 mod x^t, g(0) = g0         */ 
/*------------------------------------------------------------*/
void ZZ_pXY::series_solution(ZZ_pX & g, const ZZ_p & g0, long t){
  clear(g);
  SetCoeff(g, 0, g0);
  long prec = 1;
  ZZ_pXY dyF;
  diffY(dyF, *this);
  while (prec < t){
    prec = 2*prec;
    ZZ_pX valF, valdF, i_dF;
    eval(valF, g, prec);
    dyF.eval(valdF, g, prec);
    i_dF = InvTrunc(valdF, prec);
    valF = trunc(valF * i_dF, prec);
    g = g - valF;
  }
  g = trunc(g, t);
}


/*------------------------------------------------------------*/
/* resizes the array of coefficients                          */
/* to remove the trailing entries that are zero, if any       */
/*------------------------------------------------------------*/
void ZZ_pXY::normalize(){
  long idx = coeffX.length()-1;
  while (idx >= 0 && coeffX[idx] == 0)
    idx--;
  coeffX.SetLength(idx+1);
}

/*------------------------------------------------------------*/
/* builds from a file                                         */
/* format: [[a0,a1,..,aN],...,[x0,x1,..,xM]]                  */
/* where [a0,a1,..,aN] = coeff(f,Y^0) \in Z/pZ[X], etc.       */
/*------------------------------------------------------------*/
void ZZ_pXY::read_from_file(const string& filename){
  ifstream input;
  input.open(filename);
  string line;
  getline(input, line);
  replace(line.begin(), line.end(), ',', ' ');
  stringstream line_stream;
  line_stream.str(line);

  line_stream >> coeffX;
  input.close();
}


/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
istream& operator>>(istream& s, ZZ_pXY& x){
  NTL_INPUT_CHECK_RET(s, s >> x.coeffX);
  x.normalize();
  return s;
}

/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
ostream& operator<<(ostream& s, const ZZ_pXY& a){
  return s << a.coeffX;
}

/*------------------------------------------------------------*/
/* initializes M<y,x>=ZZ_p[y,x] with lex order y > x          */
/*------------------------------------------------------------*/
void magma_init_bi_ZZ_p(){
  cout << "M<YYY,XXX> := PolynomialRing(GF(" << ZZ_p::modulus() << ",2);\n";
}


/*------------------------------------------------------------*/
/* prints a poly in M = QQ[y > x]                             */
/*------------------------------------------------------------*/
void magma_output(const ZZ_pXY & a){
  cout << "(M!(0)";
  long d = a.degY();
  for (long i = 0; i <= d; i++){
    cout << "+(";
    magma_output(a.coeffX[i], "XXX");
    cout << ")*YYY^" << i;
   }
  cout << ")";
}

/*------------------------------------------------------------*/
/* assigns a poly in M = QQ[y > x]                            */
/*------------------------------------------------------------*/
void magma_assign(const ZZ_pXY & a, const string & name){
  cout << name << " := ";
  magma_output(a);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* derivative in Y                                            */
/*------------------------------------------------------------*/
void diffY(ZZ_pXY & dyF, const ZZ_pXY & F){
  long dy = F.degY();
  Vec<ZZ_pX> coeffs_dyF;

  if (dy <= 0){
    coeffs_dyF.SetLength(0);
  }
  else{
    coeffs_dyF.SetLength(dy);
    for (long i = 0; i < dy; i++)
      coeffs_dyF[i] = (i+1)*F.coeffX[i+1];
  }
  dyF = ZZ_pXY(coeffs_dyF);
}

/*------------------------------------------------------------*/
/* random poly with deg(F,x) <= dx, deg(F,y) <= dy            */
/*------------------------------------------------------------*/
void random(ZZ_pXY & F, long dx, long dy){
  if (dx == -1 || dy == -1){
    F.coeffX.SetLength(0);
    return;
  }

  F.coeffX.SetLength(dy+1);
  for (long i = 0; i <= dy; i++)
    random(F.coeffX[i], dx+1);
  F.normalize();
}

/*------------------------------------------------------------*/
/* reduces modulo p                                           */
/*------------------------------------------------------------*/
void conv(ZZ_pXY & f, const ZZXY & fZ){
  Vec<ZZ_pX> coeffs;
  conv(coeffs, fZ.coeffX);
  f = ZZ_pXY(coeffs);
}

#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "ZZ_pX_extra.h"
#include "ZZ_pXY.h"
#include "ZZX_extra.h"
#include "ZZXY.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* evaluate F(x,g(x)) mod x^(d^2)                             */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){
  long d = 7;
  long b = 2;

  ZZXY F;
  random(F, b, d, d);

  ZZX F0;
  random(F0, b, F.degY()-1);
  ZZX factor;
  ZZ root{22};
  SetCoeff(factor, 1, 1);
  SetCoeff(factor, 0, -root);
  F0 = F0 * factor;

  for (long i = 0; i <= F.degY(); i++)
    F.coeffX[i] = (F.coeffX[i] << 1) + coeff(F0, i);

  ZZ_pXY Fp;
  conv(Fp, F);
  
  ZZ_pX g, h;

  Fp.series_solution(g, to_ZZ_p(root), d*d);
  Fp.eval(h, g, d*d);
  cout << h << endl;

  
  // //F.eval(h, g, d*d);

  // magma_init_bi_QQ();
  // magma_init_QQX();
  // magma_assign(F, "F");
  // magma_assign(g, "XX", "g");
  // cout << "gr:=1/" << den_g << "*g;\n";
  // magma_assign(h, "XX", "h");
  // cout << "hr:=1/" << den_h << "*h;\n";
  // cout << "(hr-Evaluate(F, [gr,XX])) mod XX^" << d*d << ";\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
    
  ZZ p{9001};
  p = power(p, 10);
  ZZ_p::init(p);
    
  check(opt);
  return 0;
}

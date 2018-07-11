#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <assert.h>

#include "vec_ZZ_p_extra.h"
#include "mat_ZZ_p_extra.h"
#include "ZZ_p_toeplitz.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a toeplitz;                                    */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){
  ZZ p{9001};
  p = power(p, 4);
  ZZ_p::init(p);
  
  for (long m = 1; m < 100; m++)
    for (long n = 1; n < 100; n++){
      Vec<ZZ_p> vec;
      random(vec, m + n - 1);
      ZZ_p_toeplitz T(vec, m, n);
      Mat<ZZ_p> M;
      T.to_dense(M);
      
      Vec<ZZ_p> in, out, check;
      random(in, n);
      T.mul_right(out, in);
      mul(check, M, in);
      cout << m << " " << n << " " << (out == check) << endl;

      Mat<ZZ_p> in_m, out_m, check_m;
      long alpha = 5;
      random(in_m, n, alpha);
      T.mul_right(out_m, in_m);
      mul(check_m, M, in_m);
      cout << m << " " << n << " " << (out_m == check_m) << endl;
    }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}

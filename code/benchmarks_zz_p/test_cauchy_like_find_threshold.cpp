#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "lzz_p_cauchy_geometric.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* if opt = 1, runs FFT0                                      */
/* else, runs mod 65537                                       */
/*------------------------------------------------------------*/
void check(int opt){
  if (opt == 1)
    zz_p::FFTInit(0);
  else
    zz_p::UserFFTInit(65537);

  cout << "#p=" << zz_p::modulus() << endl;

  long step = 1;
  for (long alpha = 1; alpha < 200; alpha += step){
    if (alpha > 10)
      step = 2;
    if (alpha > 100)
      step = 5;
    long i;
    for (i = max(100, 2*alpha); ;i += 100){    
      zz_p a = to_zz_p(9);
      long j = i;
      
      mat_zz_p A, B;
      random(A, i, alpha);
      random(B, j, alpha);
      lzz_p_cauchy_like_geometric M(A, B, to_zz_p(1), power(a, i), a);
       
      long thresh = i/2-10;
      double t1, t2;
      t1 = GetTime();
      lzz_p_cauchy_like_geometric Minv;
      invert_fast(Minv, M, thresh);
      t1 = GetTime() - t1;
      
      t2 = GetTime();
      invert(Minv, M);
      t2 = GetTime() - t2;

      if (t1 < t2)
  	break;
    }
    cout << alpha << " " << i << endl;
  }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);

  return 0;
}

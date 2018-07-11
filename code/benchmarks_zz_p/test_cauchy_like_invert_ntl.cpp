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

  long lo = 100;
  long hi = 5000;

  cout << "# p=" << zz_p::modulus() << endl;

  for (long i = lo; i < hi; i += 100){

     long NB = 1;
     if (i < 1000)
       NB = 5;
     if (i < 800)
       NB = 10;
     if (i < 700)
       NB = 20;
     if (i < 600)
       NB = 60;
     if (i < 500)
       NB = 125;
     if (i < 400)
       NB = 250;
     if (i < 300)
       NB = 500;
     if (i < 200)
       NB = 1000;
     if (i < 100)
       NB = 2000;
     if (i < 50)
       NB = 5000;

     Mat<zz_p> check, check_i;
     random(check, i, i);
     double t0 = GetTime();
     for (long j = 0; j < NB; j++)
       check_i = inv(check);
     t0 = GetTime() - t0;
     t0 = t0/NB;
     cout << i << " " << 0 << " " << t0 << endl;
   }

}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if the 1st argument 1 is given, FFT0; else 65537           */
/*------------------------------------------------------------*/
int main(int argc, char** argv){
  if (argc > 1)
    check(atoi(argv[1]));

  return 0;
}

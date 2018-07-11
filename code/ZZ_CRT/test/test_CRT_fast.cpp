#include <NTL/ZZ.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "ZZ_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* does a multimod                                            */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){
  long NB = 7625;
  Vec<long> moduli;
  Vec<long> rems;
  rems.SetLength(NB);

  
  long base = 1L << 59;
  moduli.SetLength(NB);
  for (long i = 0; i < moduli.length(); i++){
    base = NextPrime(base+1);
    moduli[i] = base;
    rems[i] = rand() % base;
  }
    
  ZZ_CRT_crt_fast crt(moduli);

  ZZ b;
  double t = GetTime();
  long avg = 100;
  for (long i = 0; i < avg; i++)
    crt.crt(b, rems);
  cout << (GetTime()-t)/avg*10000 << endl;
 
} 

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}

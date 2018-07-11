#include <NTL/ZZ.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "ZZ_CRT.h"
#include "ratrecon.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long sz = 10;

  bool done = false;
  ZZ m = to_ZZ(100);
  ZZ den{1};

  while(done == false && den != 2000){
    m = NextPrime(m*10);
    Vec<ZZ> v;
    v.SetLength(3*sz);
    for (long i = 0; i < sz; i++)
      v[i] = InvMod(to_ZZ(200) % m, m);
    for (long i = sz; i < 2*sz; i++)
      v[i] = to_ZZ(0);
    for (long i = 2*sz; i < 3*sz; i++)
      v[i] = InvMod(to_ZZ(2000) % m, m);

    Vec<ZZ> w;
    long s = ReconstructRational_plain(w, den, v, m);
    if (s == 0){
      cout << "false with m=" << m << endl;
      done = false;
    }
    else{
      cout << "true with m=" << m << endl;
      done = true;
      cout << den << endl;
      cout << w << endl;
    }
  }
} 

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}

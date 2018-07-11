#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <assert.h>

#include "vec_ZZ_p_extra.h"
#include "mat_ZZ_p_extra.h"
#include "ZZ_p_block_sylvester.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* p = (2^59+1)^k                                             */
/* type = [m, ..., m] (n times)                               */
/*------------------------------------------------------------*/
void check(int opt, long k, long n, long m){

  ZZ p = to_ZZ(576460752303423489);
  p = power(p, k);
  ZZ_p::init(p);

  long prec = (m+1)*n;

  Vec<ZZ_pX> fs;
  fs.SetLength(n);
  for (long i = 0; i < n; i++)
    fs[i] = random_ZZ_pX(prec);

  Vec<long> type;
  type.SetLength(n);
  for (long i = 0; i < n; i++)
    type[i] = n;

  ZZ_p_block_sylvester_general B(fs, type, prec);

  Vec<ZZ_p> in, out, out2;
  Mat<ZZ_p> in_mat, out_mat, out_mat2;
  double t;
  cout << m << " " << n << " ";

  // check the vector version
  random(in, prec);

  t = GetTime();
  out = B.mul_left(in);
  cout << GetTime()-t << " " << n*(GetTime()-t) << " ";

  // check the matrix version
  random(in_mat, prec, n);
  t = GetTime();
  out_mat = B.mul_left(in_mat);
  cout << GetTime()-t << " ";

  if (opt == 1){
    Mat<ZZ_p> dense;
    B.to_dense(dense);
    t = GetTime();
    mul(out2, in, dense);
    cout << "(" << GetTime()-t << ") ";
    assert (out2 == out);
    t = GetTime();
    mul(out_mat2, transpose(dense), in_mat);
    assert (out_mat2 == out_mat);
    cout << "(" << GetTime()-t << ") ";
  }

  cout << endl;
}

void check(int opt){
  for (long i = 1; i < 200; i++)
    check(opt, 1, i, i);
}

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}

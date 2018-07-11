#include <NTL/vec_lzz_p.h>
#include <NTL/mat_lzz_p.h>
#include <assert.h>

#include "lzz_pX_mosaic_hankel.h"
#include "lzz_pX_CRT.h"
#include "lzz_p_cauchy_geometric.h"

NTL_CLIENT

void to_dense(Mat<zz_p>& M, const Vec<zz_p>& e){
  long n = e.length();
  M.SetDims(n, n);
  for (long i = 0; i < n; i++)
    for (long j = 0; j < n; j++)
      M[i][j] = 0;
  for (long i = 0; i < n; i++)
    M[i][i] = e[i];
}

void do_Z(Mat<zz_p>& M, long n, const zz_p& c){
  M.SetDims(n, n);
  for (long i = 0; i < n; i++)
    for (long j = 0; j < n; j++)
      M[i][j] = 0;
  for (long i = 0; i < n-1; i++)
    M[i+1][i] = 1;
  M[0][n-1] = c;
}

void random(Vec<zz_p>& v, long n, long m){
  v.SetLength(n+m-1);
  for (long i = 0; i < n+m-1; i++)
    v[i] = random_zz_p();
}

/*------------------------------------------------------------*/
/* if opt = 1, runs FFT0                                      */
/* else, runs mod 65537                                       */
/*------------------------------------------------------------*/
void check(int opt, long nblocks){
  long st, end;
  if (opt == 1){
    zz_p::FFTInit(0);
    st = 100;
    end = 5000;
  }
  else{
    zz_p::UserFFTInit(65537);
    st = 100;
    end = 5000;
  }

  cout << "# p=" << zz_p::modulus() << endl;

  for (long i = st; i < end; i += 100){

    long sz = i / nblocks;
    Vec<hankel> row;
    for (long k = 0; k < nblocks; k++){
      Vec<zz_p> dat;
      random(dat, i, sz);
      hankel h(dat, i, sz);
      row.append(h);
    }
    Vec<zz_p> dat1;
    random(dat1, i, 1);
    hankel h1(dat1, i, 1);
    row.append(h1);

    Vec< Vec<hankel> > H;
    H.SetLength(1);
    H[0] = row;
    mosaic_hankel MH(H);
    
    long NB = 10;
    // if (2*i < 1000)
    //   NB = 5;
    // if (2*i < 800)
    //   NB = 10;
    // if (2*i < 700)
    //   NB = 20;
    // if (2*i < 600)
    //   NB = 60;
    // if (2*i < 500)
    //   NB = 125;
    // if (2*i < 400)
    //   NB = 250;
    // if (2*i < 300)
    //   NB = 500;
    // if (2*i < 200)
    //   NB = 1000;
    // if (2*i < 100)
    //   NB = 2000;
    // if (2*i < 50)
    //   NB = 5000;

    lzz_p_cauchy_like_geometric CL;
    Mat<zz_p> X, Y;
    Vec<zz_p> e, f;
    zz_pX_Multipoint_Geometric X_int, Y_int;
    
    double t;

    cout << i;
    
    t = GetTime();
    for (long u = 0; u < NB; u++)
      to_cauchy_grp(CL, X_int, Y_int, e, f, MH);
    t = GetTime() - t;
    t = t/NB;
    cout << " " << t;

    t = GetTime();
    for (long u = 0; u < NB; u++){
      lzz_p_cauchy_like_geometric invCL;
      invert_fast(invCL, CL);
    }
    t = GetTime() - t;
    t = t/NB;
    cout << " " << t;

    zz_pX c, g, h, w, v;
    c = random_zz_pX(sz);
    g = random_zz_pX(sz);
    t = GetTime();
    for (long u = 0; u < NB; u++)
      XGCD(h, w, v, c, g);
    t = GetTime() - t;
    t = t/NB;
    cout << " " << t;

    cout << endl;

  }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
int main(int argc, char** argv){
  check(atoi(argv[1]), atoi(argv[2]));
  return 0;
}

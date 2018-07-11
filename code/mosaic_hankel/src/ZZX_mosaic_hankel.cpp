#include <NTL/ZZX.h>
#include <NTL/matrix.h>

#include "ZZX_mosaic_hankel.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Mosaic Hankel matrices:                            */
/* block matrix where each block is Hankel            */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void to_dense(Mat<ZZ>& Mdense, const ZZ_mosaic_hankel& M){
  long n = M.NumRows();
  long m = M.NumCols();
  Mdense.SetDims(n, m);

  for (long i = 0; i < n; i++)
    for (long j = 0; j < m; j++)
      Mdense[i][j] = M(i,j);
}

/*----------------------------------------------------*/
/* access to particular rows and columns              */
/*----------------------------------------------------*/
void first_column_of_block(Vec<ZZ>& res, long i, const ZZ_mosaic_hankel& M){
  long n = M.NumRows();
  long nb = M.NumBlockRows();

  res.SetLength(n);

  long ind = 0;
  for (long r = 0; r < nb; r++){
    const ZZ* dat = M.data[r][i].data_rev.elts();
    for (long j = 0; j < M.NumRows_of_block(r); j++)
      res[ind++] = dat[j];
  }
}

/*----------------------------------------------------*/
/* access to particular rows and columns              */
/*----------------------------------------------------*/
void last_column_of_block(Vec<ZZ>& res, long i, const ZZ_mosaic_hankel& M){
  long n = M.NumRows();
  long nb = M.NumBlockRows();

  res.SetLength(n);

  long ind = 0;
  for (long r = 0; r < nb; r++){
    const ZZ* dat = M.data[r][i].data_rev.elts();
    long shift = M.NumCols_of_block(i)-1;
    for (long j = 0; j < M.NumRows_of_block(r); j++)
      res[ind++] = dat[j + shift];
  }
}

/*----------------------------------------------------*/
/* access to particular rows and columns              */
/*----------------------------------------------------*/
void first_row_of_block(Vec<ZZ>& res, long i, const ZZ_mosaic_hankel& M){
  long m = M.NumCols();
  long mb = M.NumBlockCols();

  res.SetLength(m);

  long ind = 0;
  for (long r = 0; r < mb; r++){
    const ZZ* dat = M.data[i][r].data_rev.elts();
    for (long j = 0; j < M.NumCols_of_block(r); j++)
      res[ind++] = dat[j];
  }
}

/*----------------------------------------------------*/
/* access to particular rows and columns              */
/*----------------------------------------------------*/
void last_row_of_block(Vec<ZZ>& res, long i, const ZZ_mosaic_hankel& M){
  long m = M.NumCols();
  long mb = M.NumBlockCols();

  res.SetLength(m);

  long ind = 0;
  for (long r = 0; r < mb; r++){
    const ZZ* dat = M.data[i][r].data_rev.elts();
    long shift = M.NumRows_of_block(i)-1;
    for (long j = 0; j < M.NumCols_of_block(r); j++)
      res[ind++] = dat[j+shift];
  }
}

/*----------------------------------------------------*/
/* G, H such that Z1 M - Z0^t M = G H^t               */
/*----------------------------------------------------*/
void generators(Mat<ZZ>& G, Mat<ZZ>& H, const ZZ_mosaic_hankel& M){
  long n = M.NumRows();
  long m = M.NumCols();
  long alpha_r = M.NumBlockRows();
  long alpha_c = M.NumBlockCols();
  long alpha = alpha_r + alpha_c;

  G.SetDims(n, alpha);
  H.SetDims(m, alpha);

  Vec<ZZ> tmp_row, tmp_col;
  tmp_row.SetLength(m);
  tmp_col.SetLength(n);

  long idx;

  idx = 0;
  for (long i = 0; i < alpha_c; i++){
    first_column_of_block(tmp_col, i, M);
    G[0][i] = tmp_col[n-1];
    for (long j = 1; j < n; j++)
      G[j][i] = tmp_col[j-1];

    if (i > 0){
      last_column_of_block(tmp_col, i-1, M);
      for (long j = 0; j < n; j++)
    	G[j][i] = G[j][i] - tmp_col[j];
    }
  }
  for (long i = 0; i < alpha_r; i++){
    for (long j = 0; j < n; j++)
      G[j][i+alpha_c] = 0;
    G[idx][i+alpha_c] = 1;
    idx += M.NumRows_of_block(i);
  }


  idx = 0;
  for (long i = 0; i < alpha_c; i++){
    for (long j = 0; j < m; j++)
      H[j][i] = 0;
    H[idx][i] = 1;
    idx += M.NumCols_of_block(i);
  }
  for (long i = 0; i < alpha_r; i++){
    if (i == 0)
      last_row_of_block(tmp_row, alpha_r-1, M);
    else
      last_row_of_block(tmp_row, i-1, M);
    for (long j = 1; j < m; j++)
      H[j][i+alpha_c] = tmp_row[j];
    first_row_of_block(tmp_row, i, M);
    for (long j = 1; j < m; j++)
      H[j][i+alpha_c] = H[j][i+alpha_c] - tmp_row[j-1];
    long jdx = 0;
    for (long j = 0; j < alpha_c; j++){
      H[jdx][i+alpha_c] = 0;
      jdx += M.NumCols_of_block(j);
    }
  }
}


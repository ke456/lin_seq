#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/lzz_pEX.h>

#include "bivariate_mulmod.h"
#include "sage_output.h"

NTL_CLIENT

// checks univariate transposed multiplication
void check_transposed_product(long d){
    zz_pX a, b, c, h, inv_rev_h;
    Vec<zz_p> v;

    h = random_zz_pX(d);
    SetCoeff(h, d, 1);
    
    inv_rev_h = InvTrunc(reverse(h), d);

    a = random_zz_pX(d);
    b = random_zz_pX(d);
    
    v.SetLength(d);
    for (long i = 0; i < d; i++)
    {
	v[i] = random_zz_p();
    }
    
    c = (a * b) % h;
    zz_p res1 = to_zz_p(0);        // res1 = dotproduct(v, ab %h)
    for (long i = 0; i < d; i++)
    {
	res1 += v[i]*coeff(c, i);
    }
    
    zz_pX polV, tmp1, tmp2, tmp3;
    conv(polV, v);
    mul_Yp(tmp1, polV, h);
    tmp2 = (a * tmp1) % h;
    mul_inv_Yp(tmp3, tmp2, h, inv_rev_h);

    zz_p res2 = to_zz_p(0);     // res2 = dotproduct(a^t v, b);
    for (long i = 0; i < d; i++)
    {
	res2 += coeff(tmp3, i)*coeff(b, i);
    }
    
    cout << (res1 == res2) << endl;
}

// checks bivariate modular multiplication
// outputs results in sage format
void check_bivariate_product(long d, long m)
{
    zz_pX h;
    Vec<zz_pX> A, B, C;

    h = random_zz_pX(d);
    SetCoeff(h, d, 1);

    random(A, d, m);
    random(B, d, m);
    bivariate_product_kronecker(C, A, B, h);

    init_sage();
    sage_assign(h, "h");
    sage_assign(A, "A");
    sage_assign(B, "B");
    sage_assign(C, "C");
    cout << "I = ideal([h, x^" << m << "])\n";
    cout << "G = I.groebner_basis()\n";
    cout << "print ((C-A*B).reduce(G) == 0)\n";
}

// checks bivariate modular multiplication
// outputs results in sage format
void check_bivariate_transposed_product(long d, long m)
{
    zz_pX h;
    Vec<zz_pX> A, B, C, ell, Aell;

    h = random_zz_pX(d);
    SetCoeff(h, d, 1);

    random(A, d, m);
    random(B, d, m);
    random(ell, d, m);

    bivariate_product(C, A, B, h);
    bivariate_transposed_product(Aell, A, ell, h);
    
    zz_p res1 = to_zz_p(0);
    for (long i = 0; i < m; i++)
    {
	for (long j = 0; j < d; j++)
	{
	    res1 += coeff(C[i], j)*coeff(ell[i], j);
	}
    }

    zz_p res2 = to_zz_p(0);
    for (long i = 0; i < m; i++)
    {
	for (long j = 0; j < d; j++)
	{
	    res2 += coeff(B[i], j)*coeff(Aell[i], j);
	}
    }

    cout << (res1 == res2) << endl;
}



int main(int argc, char **argv)
{
    zz_p::init(10000000019);
    long m, d; 
    d = 4;
    m = 5;

    long test = 2;

    switch(test)
    {
    case 0:
	check_transposed_product(d);
	break;
    case 1:
	check_bivariate_product(d, m);
	break;
    case 2:
	check_bivariate_transposed_product(d, m);
	break;
    }

    return 0;
}

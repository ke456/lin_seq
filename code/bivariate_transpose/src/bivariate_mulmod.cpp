#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/vector.h>
#include <assert.h>

#include "bivariate_mulmod.h"

NTL_CLIENT

// initializes polynomials in sage
void init_sage()
{
    cout << "p = " << zz_p::modulus() << endl;
    cout << "k = GF(p)\n";
    cout << "M.<x,y> = PolynomialRing(k, order='lex')\n";
}

// prints a univariate polynomial in y
void print(const zz_pX& A)
{
    if (A == 0)
    {
	cout << "0";
	return;
    }
    for (long i = 0; i <= deg(A); i++)
    {
	cout << coeff(A, i) << "*y^" << i;
	if (i < deg(A))
	{
	    cout << " + ";
	}
    }
}

// prints a polynomial in x, with coefficients in y
void print(const Vec<zz_pX>& A)
{
    if (A.length() == 0)
    {
	cout << "0";
	return;
    }
    for (long i = 0; i < A.length(); i++)
    {
	cout << "(";
	print(A[i]);
	cout << ") * x^" << i;
	if (i < A.length()-1)
	{
	    cout << " + ";
	}
    }
}

// // prints the assignment of variable name to polynomial A (in sage)
// void sage_assign(const zz_pX& A, const string& name)
// {
//     cout << name << " = ";
//     print(A);
//     cout << endl;
// }

// prints the assignment of variable name to bivariate polynomial A (in sage)
void sage_assign(const Vec<zz_pX>& A, const string& name)
{
    cout << name << " = ";
    print(A);
    cout << endl;
}

// random vector of length m with coefficients of degree < d
void random(Vec<zz_pX>& res, long d, long m)
{
    res.SetLength(m);
    for (long i = 0; i < m; i++)
    {
	res[i] = random_zz_pX(d);
    }
}

// multiplies by the matrix Yp = [h1 h2 ... hm]
//                               [h2 ... hm 0 ]
//                               [...         ]
//                               [hm          ]
void mul_Yp(zz_pX& res, const zz_pX& a, const zz_pX& h)
{
    long m = deg(h);
    zz_pX revH = reverse(h);
    zz_pX tmp = MulTrunc(a, revH, m);
    res = reverse(tmp, m-1);
}

// multiplies by the inverse of Yp = [h1 h2 ... hm]
//                                   [h2 ... hm 0 ]
//                                   [...         ]
//                                   [hm          ]
// assumes we are given 1/rev(h) mod x^m
void mul_inv_Yp(zz_pX& res, const zz_pX& a, const zz_pX& h, const zz_pX& inv_rev_h)
{
    long m = deg(h);
    zz_pX revA = reverse(a, m-1);
    res = MulTrunc(revA, inv_rev_h, m);
}


// every body must be reduced w.r.t. h
void bivariate_product_kronecker(Vec<zz_pX>& c, const Vec<zz_pX>& a, const Vec<zz_pX>& b, const zz_pX& h)
{
    long k = deg(h);
    long dX = 2*(k-1);
    long m = a.length();
    assert (b.length() == m);
    long idx;

    zz_pX apX, bpX, cpX;
    apX.rep.SetLength((dX+1)*m);
    idx = 0;
    for (long i = 0; i < m; i++)
	for (long j = 0; j <= dX; j++)
	    SetCoeff(apX, idx++, coeff(a[i], j));
    apX.normalize();
    
    bpX.rep.SetLength((dX+1)*m);
    idx = 0;
    for (long i = 0; i < m; i++)
	for (long j = 0; j <= dX; j++)
	    SetCoeff(bpX, idx++, coeff(b[i], j));
    bpX.normalize();
    
    cpX = apX * bpX;
    
    idx = 0;
    c.SetLength(m);
    for (long i = 0; i < m; i++)
    {
	c[i] = 0;
	for (long j = 0; j <= dX; j++)
	    SetCoeff(c[i], j, coeff(cpX, idx++));
	c[i].normalize();
	c[i] %= h;
    }
}

// a and b have length EXACTLY m
// multiplies a by b modulo h(y), x^m
// result has length m
void bivariate_product(Vec<zz_pX>& res, const Vec<zz_pX>& a, const Vec<zz_pX>& b, const zz_pX& h)
{
    zz_pEBak bak;
    bak.save();
	
    long m = a.length();
    assert (b.length() == m);

    zz_pE::init(h);
    Vec<zz_pE> vecA, vecB;
    conv(vecA, a);
    conv(vecB, b);
    zz_pEX Apol, Bpol, Rpol;
    conv(Apol, vecA);
    conv(Bpol, vecB);

    Rpol = MulTrunc(Apol, Bpol, m);
    
    res.SetLength(m);
    for (long i = 0; i < m; i++)
	conv(res[i], coeff(Rpol, i));

    bak.restore();
}

// a and b have length EXACTLY m
// res = a^t b
// result has length m
void bivariate_transposed_product(Vec<zz_pX>& res, const Vec<zz_pX>& a, const Vec<zz_pX>& b, const zz_pX& h){
    zz_pX inv_rev_h;
    Vec<zz_pX> revB, tmp;
    long m = a.length();
    assert (b.length() == m);
    long d = deg(h);

    inv_rev_h = InvTrunc(reverse(h), d);

    revB.SetLength(m);
    for (long i = 0; i < m; i++)
    {
	mul_Yp(revB[i], b[m-1-i], h);
    }

    bivariate_product_kronecker(tmp, a, revB, h);

    res.SetLength(m);
    for (long i = 0; i < m; i++)
    {
	mul_inv_Yp(res[i], tmp[m-1-i], h, inv_rev_h);
    }
}

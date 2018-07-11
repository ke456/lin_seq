#include "lin_seq.h"
#include "bivariate_mulmod.h"
#include <iostream>
using namespace std;

void lin_seq::sqfree_decomp(vec_pair_zz_pX_long& result,const zz_pX f){
	auto f_d = diff(f);
	auto u = GCD(f,f_d);
	auto v = f/u;
	auto w = f_d/u;
	
	result.SetLength(0);

	long i  = 1;
	
	do{
		auto v_d = diff(v);
		auto h = GCD(v, w-v_d);
		if (deg(h) > 0){
			append(result, cons(h,i));
		}
		i++;
		v = v/h;
		w = (w-v_d)/h;
	}while(v != zz_pX(1));
}

void inverse_mod_newton(zz_pX &result, const zz_pX& f, zz_pX g, zz_pX P, const long m){
	auto r = ceil(long(m)/long(2));
	for (long i = 1; i <= r; i++){
		mul(P,P,P);
		g = (2*g-f*g*g)%P;
	}
	result = g;
}

void lin_seq::hankel_solve(zz_pX &result, const Vec<zz_p> &lhs, 
                           const Vec<zz_p> &rhs, const long ind){
	auto &P = expanded[ind];
	zz_pX n_lhs, n_rhs;
	//cout << "lhs: " << lhs << endl;
	//cout << "rhs: " << rhs << endl;
	auto n = rhs.length();
	for (long i = 0; i<n; i++){
		SetCoeff(n_lhs, i, lhs[i]);
		SetCoeff(n_rhs, i, rhs[i]);
	}
	
	auto denom = reverse(P);	

	// compute the numerators
	MulTrunc(n_lhs,n_lhs,denom,n);
	MulTrunc(n_rhs,n_rhs,denom,n);
	//cout << "n_lhs: " << n_lhs << endl;
	//cout << "n_rhs: " << n_rhs << endl;

	// reverse each polynomial
	reverse(n_lhs,n_lhs, n-1);
	reverse(n_rhs,n_rhs, n-1);

	zz_pX inv;
	InvMod(inv, n_lhs, P);
	MulMod(result, n_rhs, inv, P);

}

lin_seq::lin_seq(const zz_pX& P, const Vec<zz_p>& init, const long &p):
	p{p}, P{P}, init{init}{
	auto t = GetTime();
	
	sqfree_decomp(sq_free_pols,P);
	cout << "sq_free 2 decomp took: " << GetTime() - t << endl;
	cout << "sq_free decom: " << sq_free_pols << endl;

	expanded.SetLength(sq_free_pols.length());
	sq_free_mods.SetLength(sq_free_pols.length());
	t = GetTime();
	long max_md = 0;
	for (long i = 0; i < sq_free_pols.length(); i++){
		repeated_sq(expanded[i], sq_free_pols[i].a, ZZ(sq_free_pols[i].b));
		if (sq_free_pols[i].b*deg(sq_free_pols[i].a) > max_md)
			max_md = sq_free_pols[i].b*deg(sq_free_pols[i].a);
		build(sq_free_mods[i],sq_free_pols[i].a);
	}
	cout << "expanding took: " << GetTime() - t << endl;
	crt = new zz_pX_CRT(expanded);
	
	// precompute the factorials up to max_md
	//	cout << "max_md: " << max_md << endl;
	factorials.SetLength(max_md+1);
	inv_factorials.SetLength(max_md+1);
	factorials[0] = inv_factorials[0] = zz_p{1};
	for (long i = 1; i <= max_md; i++){
		factorials[i] = factorials[i-1]*i;
		//		cout << i<< "!: " << factorials[i] << endl;
		inv_factorials[i] = 1/factorials[i];
	}
	
}

void lin_seq::gen_binomial(vec_zz_p& v, const long m, const ZZ &D){
	v.SetLength(m);
	v[0] = 1;
	zz_p d(rem(D,p));
	for (long i = 0; i < m; i++)
		v[i+1] = (d-i) / (i+1) * v[i];
}

void lin_seq::repeated_sq(zz_pX &result, const zz_pX &p, const ZZ &n){
	zz_pX t = zz_pX(1);
	auto x = p;
	ZZ k = n;
	while (k > ZZ(0)){
		if (k%ZZ(2) == 1){
			mul(t,t,x);
			k -= ZZ(1);
		}else{
			mul(x,x,x);
			k /= ZZ(2);
		}
	}
	result = t;
/*
	cout <<"sq1: " << result << endl;


	auto bin = get_binary(n);
	zz_pX b = p;
	for (long i = bin.length() - 2; i >= 0; i--){
		mul(b,b,b);
		if (bin[i] == 1)
			mul(b,b,p);
	}
	result = b;
	cout <<"sq2: " << result << endl;
*/
}

void lin_seq::repeated_sq_mod(zz_pX &result, const zz_pX &p,
                              const zz_pX &h, const ZZ &n){
	zz_pX t = zz_pX(1);
	zz_pX x = p%h;
	zz_pX b = p % h;
	zz_pX p_mod = b;

	ZZ k = n;
	while (k > ZZ(0)){
		if (k%ZZ(2) == 1){
			MulMod(t,t,x,h);
			k -= ZZ(1);
		}else{
			MulMod(x,x,x,h);
			k /= ZZ(2);
		}
	}
	rem(t,t,h);
	result = t;

/*	
	auto bin = get_binary(n);
	for (long i = bin.length()-2; i>=0; i--){
		MulMod(b,b,b,h);
		if (bin[i] == 1)
			MulMod(b,b,p_mod,h);
	}
	result = b;
*/
}

Vec<int> lin_seq::get_binary(const ZZ &t){
  Vec<int> result;
	auto n = t;
	while (n != 0){
		result.append(n % 2);
		n = n/2;
	}
	return result;
}

void lin_seq::compute_x(zz_pX &result, const ZZ &D, const long m){
	auto t = GetTime();
	vec_zz_p bin_coeffs;
	gen_binomial(bin_coeffs,m,D);
	
	// result = (x+1)^D % x^m
	result = zz_pX();
	for (long i = 0; i < m; i++)
		SetCoeff(result,i,bin_coeffs[i]);
	cout << "bin coeff took: " << GetTime()-t << endl;
	t = GetTime();
	shift_by(result, result, zz_p(-1));
	cout << "mul shift took: " << GetTime() -t  << endl;
}

zz_p fac (long n){
	zz_p result{1};

	for (long i = 1; i <= n; i++)
		result *= i;
	
	return result;
}

Vec<zz_p> facs (long n){
	Vec<zz_p> result;
	result.SetLength(n+1);
	result[0] = zz_p{1};
	for (long i = 1; i <= n; i++)
		result[i] = result[i-1] *i;
	return result; 
}

void lin_seq::get_h_p(zz_pX &result, const long &n){
	result = zz_pX();
	for (long i = 0; i < n; i++)
		SetCoeff(result, i, inv_factorials[i]);	
}

void lin_seq::mult_shift(Vec<zz_p> &result, const Vec<zz_p> &v){
	auto n = v.length();
	result.SetLength(n);
	zz_pX rhs;
	for (long i = 0; i < n; i++){
		SetCoeff(rhs, n-i-1, factorials[i] * v[i]);
	}
	zz_pX hp;
	get_h_p(hp, n);
	MulTrunc(rhs, hp, rhs, n);

	for (long i = 0; i < n; i++){
		result[i]=(inv_factorials[i] * coeff(rhs,n-1-i));
	}
}

void lin_seq::T_mult_shift(Vec<zz_p> &result, const Vec<zz_p> &v){
	auto n = v.length();
	result.SetLength(n);
	zz_pX rhs;
	for (long i = 0; i < n; i++)
		SetCoeff(rhs, i, 1/fac(i) * v[i]);
	zz_pX hp;
	get_h_p(hp,n);
	MulTrunc(rhs, hp, rhs, n);
	for (long i = 0; i < n; i++)
		result[i] = (fac(i) * coeff(rhs,i));
}

void lin_seq::shift(zz_pX &result, const zz_pX &P){
	Vec<zz_p> rhs;
	for (long i = 0; i <= deg(P); i++){
		rhs.append(coeff(P,i));
	}
	mult_shift(rhs,rhs);
	result = zz_pX();
	for (long i = 0; i < rhs.length(); i++)
		SetCoeff(result, i, rhs[i]);
}

void lin_seq::T_shift(zz_pX &result, const zz_pX &P){
	Vec<zz_p> rhs;
	for (long i = 0; i <= deg(P); i++){
		rhs.append(coeff(P,i));
	}
	T_mult_shift(rhs,rhs);
	result = zz_pX();
	for (long i = 0; i < rhs.length(); i++)
		SetCoeff(result, i, rhs[i]);
}

void lin_seq::shift_by(zz_pX &result, const zz_pX &P, const zz_p &s){
	result = P;
	
	Vec<zz_p> pow;
	pow.SetLength(deg(P)+1);
	pow[0] = zz_p(1);
	Vec<zz_p> inv_pow;
	inv_pow.SetLength(deg(P)+1);
	inv_pow[0] = zz_p(1);

	for (long i = 1; i <= deg(P); i++){
		pow[i] = pow[i-1] * s;
		inv_pow[i] = 1/pow[i];
	}
	for (long i = 0; i <= deg(P); i++){
		SetCoeff(result,i,coeff(result,i)*pow[i]);
	}
	shift(result, result);
	for (long i = 0; i <= deg(result); i++)
		SetCoeff(result, i, coeff(result,i)*(inv_pow[i]));
}

void lin_seq::T_shift_by(zz_pX &result, const zz_pX &P, const zz_p &s){
	result = P;
	for (long i = 0; i <= deg(result); i++)
		SetCoeff(result, i, coeff(result,i)*(1/power(s,i)));
	shift(result, result);
	for (long i = 0; i <= deg(P); i++){
		SetCoeff(result,i,coeff(result,i)*power(s,i));
	}
}

Vec<zz_p> convert (const zz_pX &f){
	Vec<zz_p> r;
	r.SetLength(deg(f)+1);
	
	for (long i = 0; i <= deg(f); i++)
		r[i] = coeff(f,i);
	
	return r;
}

zz_pX convert (const Vec<zz_p> &v){
	zz_pX r;
	for (long i = 0; i < v.length(); i++)
		SetCoeff(r,i,v[i]);
	return r;
}

// bivariate shifts
void lin_seq::shift_by_alpha(zz_pXY &result, const zz_pXY &f, 
                             const long ind, const int dir){
	result = f;
	long m = result.coeffX.length();
	auto &h = sq_free_pols[ind].a;
	zz_pX y; // will keep track of powers of y
	y.SetLength(m);
	SetCoeff(y,0,1); // initially set to y^0
	//	cout << "h: " << h << endl;
	// multiply by alpha^i = y^i
	for (long i = 0; i < m; i++){
		MulMod(result.coeffX[i],result.coeffX[i], y, h);
		y *= dir;
		rem(y,LeftShift(y,1),h);
	}
	//	cout << "after y^i: " << result << endl;
	// multiply by diag(i!)
	for (long i = 0; i < m; i++)
		result.coeffX[i] = result.coeffX[i] * factorials[i];
	//	cout << "after i!: " << result << endl;
	// reverse f
	for (long i = 0; i < m; i++)
		reverse(result.coeffX[i],result.coeffX[i],deg(h)-1);
	auto temp = result;
	for (long i = 0; i < m; i++)
		result.coeffX[i] = temp.coeffX[m-i-1];
	//	cout << "after reverse: " << result << endl;
		
	zz_pX shift_poly;
	get_h_p(shift_poly, m);
	zz_pXY shift_bi_poly;
	shift_bi_poly.coeffX.SetLength(deg(shift_poly)+1);
	//	cout << "shift poly: " << shift_poly << endl;
	for (long i = 0; i <= deg(shift_poly); i++){
		shift_bi_poly.coeffX[i] = zz_pX(coeff(shift_poly,i));
	}
	//	cout << "bi-shift: " << shift_bi_poly << endl;
	mul_kronecker(result, result, shift_bi_poly);
	//	cout << "mult: " << result << endl;
	result.coeffX.SetLength(m);
	temp = result;
	for (long i = 0; i < m; i++)
		result.coeffX[i] = reverse(temp.coeffX[m-i-1], deg(h)-1);
	//	cout << "after trunc: " << result << endl;

	// post multiply by 1/y^i
	y = zz_pX{0};
	SetCoeff(y,1,dir);
	rem(y,y,h);
	auto inv_y = InvMod(y,h);
	SetCoeff(y,1,0);
	SetCoeff(y,0,1);
	for (long i = 0; i < m; i++){
		MulMod(result.coeffX[i],result.coeffX[i],y,h);
		//		cout << "mod h: " << result.coeffX[i] << endl;
		MulMod(y,y,inv_y,h);
		//		cout << "y: " << y << endl;
	}
	// post multiply by 1/i!
	for (long i = 0; i < m; i++)
		result.coeffX[i] = inv_factorials[i]*result.coeffX[i];
}

void lin_seq::T_shift_by_alpha(zz_pXY &result, const zz_pXY &f, const long ind, const int dir){
	result = f;
	long m = result.coeffX.length();
	auto &h = sq_free_pols[ind].a;
	auto &h_mod = sq_free_mods[ind];
	// pre-multiply by 1/i!
	for (long i = 0; i < m; i++)
		result.coeffX[i] = inv_factorials[i]*result.coeffX[i];
	
	// pre-multiply by 1/y^i transpose
	zz_pX y;
	SetCoeff(y,1,dir);
	rem(y,y,h);
	auto inv_y = InvMod(y,h);
	//	cout << "inv_y: " << inv_y << endl;
	SetCoeff(y,1,0); 
	SetCoeff(y,0,1);
	for (long i = 0; i < m; i++){
	  //		cout << "power of inv_y: " << y << endl;
		zz_pXMultiplier B(y, h);
		auto right_prod = UpdateMap(convert(result.coeffX[i]), B, h_mod);
		// cout << "vec: " << result.coeffX[i] << endl;
		// cout << "right_prod: " << right_prod << endl;
		MulMod(y, y, inv_y, h);
		result.coeffX[i] = convert(right_prod);	
	}
	
	// Toeplitz multiplication
	zz_pX shift_poly;
	get_h_p(shift_poly, m);
	zz_pXY shift_bi_poly;
	shift_bi_poly.coeffX.SetLength(deg(shift_poly)+1);
	for (long i = 0; i<=deg(shift_poly);i++)
		shift_bi_poly.coeffX[i] = zz_pX(coeff(shift_poly,i));
	// cout << "bi-shift: " << shift_bi_poly << endl;
	mul_kronecker(result, result, shift_bi_poly);
	// cout << "mult: " << result << endl;
	result.coeffX.SetLength(m);
	// cout << "after trunc: " << result << endl;
	
	// post-multiply by i!
	for (long i = 0; i < m; i++)
		result.coeffX[i] = result.coeffX[i] * factorials[i];
	// post-multiply by y^i transpose
	y = zz_pX();
	SetCoeff(y,1,dir);
	rem(y,y,h);
	zz_pX pow_y = zz_pX(zz_p(1));
	for (long i = 0; i < m; i++){
	  //		cout << "y: " << pow_y << endl;
		zz_pXMultiplier B(pow_y,h);
		auto right_prod = UpdateMap(convert(result.coeffX[i]), B, h_mod);
		//		cout << "vec: " << result.coeffX[i] << endl;
		//		cout << "right_prod: " << right_prod << endl;
		MulMod(pow_y,y,pow_y,h);
		result.coeffX[i] = convert(right_prod);	
	}
}

void lin_seq::t_expand(Vec<zz_p> &result, const zz_pX &a, const zz_pX &p, const ZZ &k){
	if (k == 1){
		result.append(coeff(a,0));
		return;
	}
	auto t = k / 2;
	zz_pX pt;
	repeated_sq(pt,p,t);

	auto q = a / pt;
	auto r = a % pt;

	t_expand(result, r, p, t);
	t_expand(result, q, p, t);
}

void lin_seq::bivariate_mod (Vec<zz_pX> &result, const ZZ &D,
                             const zz_pX &h, const long m){
	zz_pX mx;
	compute_x(mx, D, m);
	
	zz_pX y; // represents polynomial for p(y) = y
  SetCoeff(y,1,1);
	rem(y,y,h);
	zz_pX ypow;
	repeated_sq_mod(ypow, y, h, D-m+1);
	result.SetLength(m);
	for (long i = 0; i < m; i++){
		result[m-i-1] = coeff(mx,m-i-1) * ypow;
		//		cout << "result: " << result[m-i-1] << endl;
		MulMod(ypow, ypow, y, h);
	}
	//	cout << "bivariate mod: " << result << endl;
}


void untangling_hbar_powers(long m, const zz_pX &hbar, 
														map<long,zz_pX>& out){
	if(m==1){
		out[1] = hbar;
		return;
	}
	long l = m / 2;
	//	cout << "l: " << l << endl;
	if(out.find(m) == out.end()){ 
		//check if the power has already been inserted
		untangling_hbar_powers(l,hbar,out);
		untangling_hbar_powers(m-l,hbar,out);
		if((m%2)==0){
			out[m] = out[l]*out[l];
		} else {
			out[m] = out[l]*out[m-l];
		}
	}
}

void lin_seq::untangling_transpose(zz_pX &result, const zz_pXY &in,
                                   const long ind){
	zz_pXY shift_in = in;
	T_shift_by_alpha(shift_in, in, ind, -1);
	// prepare for untangling:
	auto &h = sq_free_pols[ind].a;
	auto &m = sq_free_pols[ind].b;
	map<long, zz_pX> h_pow;
	untangling_hbar_powers(m,h,h_pow);	
	for (long i = 0; i < shift_in.coeffX.length(); i++)
	  shift_in.coeffX[i] *= inv_factorials[i];


	result = untangling_transpose(shift_in, h_pow, m);
}

zz_pX lin_seq::untangling_transpose(zz_pXY u, const map<long,zz_pX>& hbar_powers, long m){
	//cout << "u: " << u << endl;
	if( m==1 ){
		//cout << "returning: " << u << endl;
		return u.coeffX[0];
	}

	auto lth_derivative_transpose = [&](long l, zz_pX &target) -> void {
		//cout << "Derivative with l=" << l << endl;
		//cout << "in: " << target << endl;
		target = LeftShift(target, l);
		for(int p=1; p<deg(target); p++){ 
			target[p] = target[p]*factorials[p]*inv_factorials[p-l];
		}
		//cout << "out: " << target << endl;
	};

	auto modulo_transpose = [&](zz_pX target, zz_pX h_power, long m) ->  zz_pX {
		target = MulTrunc(target, reverse(h_power), deg(h_power)); //numerator
		zz_pX hpower_inv = InvTrunc(reverse(h_power), m);
			return MulTrunc(target, hpower_inv, m);
	};

	long l = m / 2;
	zz_pXY first_half, second_half;
	first_half.coeffX.SetLength(l);
	second_half.coeffX.SetLength(m-l);
	for(long i=0; i<m; i++){
		if(i<l){
			first_half.coeffX[i] = u.coeffX[i];
		}else{
			second_half.coeffX[i-l] = u.coeffX[i];
		}
	}
	auto d = deg(hbar_powers.at(1));
	
	zz_pX r0 = untangling_transpose(first_half, hbar_powers, l);
	zz_pX s0 = modulo_transpose(r0, hbar_powers.at(l), m*d);
	//cout << "s0: " << s0 << endl;
	zz_pX r1 = untangling_transpose(second_half, hbar_powers, m-l);
	zz_pX s1 = modulo_transpose(r1, hbar_powers.at(m-l), m*d);
	//cout << "before derivative: " << s1 << endl;
	lth_derivative_transpose(l,s1);
	trunc(s1, s1, m*d);
	//cout << "after derivative: " << s1 << endl;
	//cout << "before return " << (s0+s1) << endl;
	return s0+s1;
}

void lin_seq::mod(zz_pX& result, const ZZ &D, const long i){
	auto &sq_free = sq_free_pols[i];
	auto &h = sq_free.a;
	auto &m = sq_free.b;

	double t = GetTime();
	Vec<zz_pX> r;
	bivariate_mod(r, D, h, m);
	cout << "bivariate mod took: " << GetTime()-t<< endl;

	t = GetTime();
	// make random linear form
	zz_pXY lambda;
	lambda.coeffX.SetLength(m);
	for (long i = 0; i < m; i++)
		for (long j = 0; j < deg(h); j++)
			SetCoeff(lambda.coeffX[i],j, random_zz_p());
	//cout << "linear form: " << lambda << endl;
	zz_pX h_seq;
	untangling_transpose(h_seq, lambda, i);
	//cout << "lhs seq: " << h_seq << endl << endl;
	cout << "lhs seq took: " << GetTime() - t << endl;
	
	t = GetTime();
	// do transpose product
	zz_pXY A_s;
	shift_by_alpha(A_s,r,i,1);
	zz_pXY L_s;
	T_shift_by_alpha(L_s,lambda, i, -1);
	zz_pXY B_s;
	bivariate_transposed_product(B_s.coeffX, A_s.coeffX, L_s.coeffX, h);
	cout << "transpose prod took: " << GetTime() - t << endl;

	t = GetTime();
	T_shift_by_alpha(B_s,B_s,i,1);
	cout << "shift took: " << GetTime() - t << endl;
	zz_pX rhs_seq;
	untangling_transpose(rhs_seq, B_s, i);
	cout << "rhs seq took: " << GetTime()-t << endl;
	//cout << "rhs seq: " << rhs_seq << endl;
	
	//zz_pX temp_result;
	t = GetTime();
	hankel_solve(result, convert(h_seq), convert(rhs_seq), i);
	cout << "hankel solve took: " << GetTime()-t << endl;
	//cout << "result: " << result << endl;
	//cout << "temp result: " << temp_result << endl;
}

void lin_seq::get_entry(zz_p& result, const ZZ &D){
	Vec<zz_pX> mods;
	mods.SetLength(sq_free_pols.length());
	for (long i = 0; i < sq_free_pols.length(); i++){
		mod(mods[i], D, i);
	}
	double t = GetTime();
	zz_pX combined;
	crt->combine(combined, mods);
	cout << "crt took: " << GetTime()-t << endl;
	result = 0;
	for (long i = 0; i < init.length(); i++){
		result = result + coeff(combined,i)*init[i];
	}
}

void lin_seq::get_entry_naive(zz_p& result, const ZZ &D){
	zz_pX mod;
	SetCoeff(mod,1,1);
	repeated_sq_mod(mod, mod, P, D);
	result = 0;
	for (long i = 0; i < init.length(); i++)
			result = result + coeff(mod,i) * init[i];
}

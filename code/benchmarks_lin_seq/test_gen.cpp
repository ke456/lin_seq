#include <iostream>
#include <fstream>
#include <sstream>
#include <NTL/lzz_pX.h>

NTL_CLIENT

using namespace std;

Vec<int> get_binary(const ZZ &t){
  Vec<int> result;
	auto n = t;
	while (n != 0){
		result.append(n % 2);
		n = n/2;
	}
	return result;
}

void repeated_sq(zz_pX &result, const zz_pX &p, const ZZ &n){
	auto bin = get_binary(n);
	zz_pX b = p;
	for (long i = bin.length() - 2; i >= 0; i--){
		mul(b,b,b);
		if (bin[i] == 1)
			mul(b,b,p);
	}
	result = b;
}

int main(){
	cout << "Enter p,d,i: ";
	long p, d, i;
	long m = 1;
	cin >> p >> d >> i;
	zz_p::init(p);
	zz_pX P{zz_p{1}};
	for(long at = 0; at < i; at++){
		m *= 2;
	}
	zz_pX f = random_zz_pX(d+1);
	repeated_sq(f,f,ZZ(m));
	mul(P,P,f);

	MakeMonic(P);
	
	ostringstream oss;
	oss << "t_"<<p<<"_"<<d<<"_"<<i;
	ofstream ofs(oss.str());
	ofs << p << endl;
	for(long at = 0; at < deg(P)+1; at++)
		ofs << coeff(P,at) << " ";
	ofs << endl;
	for (long at = 0;  at< deg(P); at++)
		ofs << 1 << " ";
	ofs << endl;
	for (long at = 1; at <= 10; at++)
		ofs << at*100*deg(P) << endl;
	
		
}







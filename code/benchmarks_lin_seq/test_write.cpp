#include "lin_seq.h"
#include <iostream>
#include <string>
#include <sstream>
#include "lzz_pXY.h"
using namespace std;

//#define hankel_test
#define lin_seq_test

int main(){
#ifdef hankel_test
	zz_p::init(13);
	cout << "STARTING HANKEL SOLVER TEST" << endl;
	long in;
	string s;
	istringstream iss;
	
	
	cout << "Enter lhs (one line):" << endl;
	Vec<zz_p> lhs;
	getline(cin,s);
	iss = istringstream(s);
	while (iss >> in){
		lhs.append(zz_p(in));
	}
	cout << "lhs: " << lhs << endl;

	cout << "Enter rhs (one line):" << endl;
	Vec<zz_p> rhs;
	getline(cin,s);
	iss = istringstream(s);
	while (iss >> in)
		rhs.append(zz_p(in));
	cout << "rhs: " << rhs << endl;

	cout << "Enter P (one line):" << endl;
	zz_pX Px;
	long i = 0;
	getline(cin,s);
	iss = istringstream{s};
	while (iss >> in)
		SetCoeff(Px, i++, zz_p(in));
	
	zz_pX result;
	lin_seq::hankel_solve(result, lhs, rhs, Px);
	cout << "result: " << result << endl;

	return 0;

#endif
#ifdef lin_seq_test
  cout << "Enter p:" << endl;
  long p;
  cin >> p;
  zz_p::init(p);
  long prime = p;

  string l;
  getline(cin,l);
  zz_pX P;
  cout << "Enter P (space separated, one line, increasing in degree)" << endl;
  getline(cin,l);
  istringstream iss(l);
  long i = 0;
  while (iss >> p){
    SetCoeff(P,i++,zz_p(p));
  }

  Vec<zz_p> init;
  cout << "Enter init conditions (space separated, one line)" << endl;
  getline(cin,l);
  iss = istringstream(l);
  while (iss >> p)
    init.append(zz_p(p));
	cout << endl;

	double t = GetTime();
  lin_seq q{P,init,prime};
	cout << "Ctor: " << GetTime() - t << endl;
  zz_p r; 
  ZZ D;
  cout << "Enter D: " << endl;
  while (cin >> D){
		cout << "D: "<<D << endl;
		t = GetTime();
    q.get_entry(r, D);
    cout << "D-th entry: " << r << endl;
		cout << "Time: " << GetTime() - t << endl <<  endl;
		t = GetTime();
		q.get_entry_naive(r,D);
		cout << "Naive: " << r << endl;
		cout << "took: " << GetTime()-t<<endl<<endl;
    cout << "Enter D: ";
  }
	cout << endl;
#endif
}











//============================================================================
// Name        : lifttest.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

double op_sum(double i, double j) {
        return (i+j);
}

int vecsum(vector<double> &a, vector<double> &b, vector<double> &c){


    c.resize(a.size());
        transform(a.begin(), a.end(), b.begin(), b.begin(), op_sum);
        c = b;
                return 0;
}

int sign(int X) {
        if (X >= 0)
                return 1;
        else
                return -1;
}

void circshift(vector<int> &sig_cir, int L){
        if ( abs(L) >(signed int) sig_cir.size()) {
                L = sign(L) * (abs(L) % sig_cir.size());
        }

        if ( L < 0 ){
                L = (sig_cir.size() + L) % sig_cir.size();
        //	cout << "L" << L << endl;
        }
                for (int i = 0; i < L; i++){
                        sig_cir.push_back(sig_cir[0]);
                        sig_cir.erase(sig_cir.begin());
                }

}

void haarlift(const vector<double> &inp) {
	vector<double> sige;
	vector<double> sigo;
	vector<double> sig = inp;

	int len=sig.size();


	if ( 2 * (int) (len/2) != len) {

		sig.push_back(0.0);
		len=len+1;
	}
   // Forward Transform

	for (int i=0; i < len/2; i++) {
		double temp1 = sig[2*i];
		double temp2 = sig[2*i+1];
		sige.push_back(temp1);
		sigo.push_back(temp2);

		sigo[i]=sigo[i]-sige[i];
		sige[i]=sige[i]+sigo[i]/2.0;
		cout << "diff" << sigo[i] << endl;
		cout << "avg" << sige[i] << endl;
	}

	// Inverse Transform
	vector<double> oup;


	for (int i=0; i < len/2; i++) {
		sige[i]=sige[i]-sigo[i]/2.0;
		sigo[i]=sige[i]+sigo[i];
		oup.push_back(sige[i]);
		oup.push_back(sigo[i]);
		cout << oup[2*i] << " " << oup[2*i+1] << " " ;
		cout << sige.size() << " " << sigo.size() << endl;
	}



}

int gcd(int a, int b) {
	/*
	 function f=euclid(a,b)
% Euclidean Algorithm
if a>b
    p(1)=a;
    p(2)=b;
else
    p(1)=b;
    p(2)=a;
end
i=1;
p(3)=1;
while p(i+1)~=0
p(i+2)=mod(p(i),p(i+1));
i=i+1;
end
'Greatest Common Divisor'
j=p(i)
	 */
	int p;
	vector<int> vec;
	if (a > b) {
		vec.push_back(a);
		vec.push_back(b);
	} else {
		vec.push_back(b);
		vec.push_back(a);
	}
	int i=0;
	int t=vec[1];
	vec.push_back(1);

	while (t != 0) {
		vec[i+2]=vec[i]%vec[i+1];
		i++;
		t=vec[i+1];
	}
	p=vec[i];
	return p;

}

int degree(vector<int> &vec) {
	int l=vec.size();
	int i=0;
	while (vec[i]==0) {
		i++;
		}

	return l-i;
}

void polydiv(const vector<int> &num, const vector<int> &den) {
	vector<int> num1=num;
	vector<int> den1=den;
	//vector<int> den2=den;

	int it=0;


	while (den1[it]==0) {
		den1.erase(den1.begin());
	}

	int degd=degree(den1);
	int lnum=num1.size();
	int lden=den1.size();
	for (int i=0;i < lnum-lden;i++) {
		den1.insert(den1.end(),0);
	}

	vector<int> quot;
	vector<int> rem;
	for (int i=0;i<lnum;i++) {
		//quot.push_back(0);
		rem.push_back(num1[i]);
		cout << den1[i] << endl;
	}
	int iter=0;
	while (degree(rem) >= degd) {
		int q1=rem[iter]/den1[iter];
		iter++;

		for (int i=0;i<lnum;i++) {
			rem[i]=rem[i]-q1*den1[i];
		}
		circshift(den1,-1);
		quot.push_back(q1);
	}

	for (int i=0;i<lnum;i++) {
				cout << rem[i] << endl;
			}

	for (unsigned int i=0;i<quot.size();i++) {
					cout << quot[i] << endl;
				}




}


int main() {
	double re_array[]= {0.8580,1.2540,-1.5937,-1.4410,0.5711,0.4700,-0.3999,0.6900};

	vector<double> a(re_array, re_array + sizeof(re_array)/sizeof(double));
	haarlift(a);
	int x=182;
	int y=546;

	int oup=gcd(x,y);
	int n1[]={4,3,0,2,1};
	int d1[]={0,0,1,1,2};
	//vector<int> num,den;
	vector<int> num(n1, n1 + sizeof(n1)/sizeof(int));
	vector<int> den(d1, d1 + sizeof(d1)/sizeof(int));
	polydiv(num,den);

	//circshift(a,1);
	//cout << a[0] << endl; // prints FFT Test
	return 0;
}

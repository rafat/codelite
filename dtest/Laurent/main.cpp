//============================================================================
// Name        : Laurent.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   :
// Description : 
//============================================================================

#include <iostream>
#include <vector>
#include "Laurent.h"

using namespace std;

int main() {
//    double re_array[]= {5.0,2.0,4.0,3.0};
//	double re_array2[]= {2.0,2.0,-1.0,9.0};
    double re_array[]= {-0.125,0.75,-0.125};
	double re_array2[]= {0.25,0.25};

	vector<double> a(re_array, re_array + sizeof(re_array)/sizeof(double));
	vector<double> b(re_array2, re_array2 + sizeof(re_array2)/sizeof(double));
//	int low=-2;
//	int low2=2;
    int low=1;
	int low2=1;
	Laurent<double> lp1,lp2,lp3,lp4,lp5,lp6;
	lp1.setPoly(a,low);
	lp2.setPoly(b,low2);
	cout << lp1.degree()+1 << " " << lp2.degree()+1 << endl; // prints Laurent Test

    vector<double> newvec;
    lp1.getPoly(newvec);
    cout << newvec.size() << endl;

    cout << "A : ";
	lp1.dispPoly();
	cout << "B : ";
	lp2.dispPoly();
	lp3.LaurentAdd(lp1,lp2);
	cout << "A+B : ";
	lp3.dispPoly();
	lp3.LaurentSub(lp1,lp2);
	cout << "A-B : ";
    lp3.dispPoly();

    lp4.LaurentMult(lp1,lp2);
    cout << "A*B : ";
    lp4.dispPoly();

    lp5.zinv(lp1);
    cout << "A(1/z) : ";
    lp5.dispPoly();

    //cout << lp5.isMono() << endl;

    LaurentMat<double> Mat1;
    Mat1.setMat(lp1,lp1,lp2,lp2);
    Laurent<double> det1;
    Mat1.Det(det1);
    cout << "Determinant : ";
    det1.dispPoly();

    cout << "Mono : " << det1.isMono() << endl;
	vector<Laurent<double> > lout;
    Div(lp1,lp2,lout);

    for (unsigned int i=0; i < lout.size()/2; i++) {
        cout << "Q" << i << ": ";
        lout[2*i].dispPoly();
        cout << endl;
        cout << "R" << i << ": " ;
        lout[2*i+1].dispPoly();
        cout << endl;

    }
	
	return 0;
}

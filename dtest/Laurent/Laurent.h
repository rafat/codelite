/*
 * Laurent.h
 *
 *  Created on: Jul 22, 2012
 *      Author: Rafat Hussain
 */

#ifndef LAURENT_H_
#define LAURENT_H_
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

template <typename T>
class Laurent {
	int deg,highest;
	vector<T> poly;

public:
	Laurent() {
	int degree=0;
	vector<T> poly;
	int highest = 0;
	// TODO Auto-generated constructor stub

}
	void setPoly(const vector<T> &coef, int hdeg) {
	deg=coef.size()-1;
	poly = coef;
	highest = hdeg;
	// TODO Auto-generated constructor stub

}
    void getPoly(vector<T> &vec) {
        vec=poly;
    }

	void dispPoly() {
	int sz=poly.size();

	for (int i = 0;i < sz; i++) {
		if (i == sz-1) {
			cout << poly[i] << "*z^" << "(" << highest-i << ")";
		}
		else if (poly[i+1] >= 0) {
			cout << poly[i] << "*z^" << "(" << highest-i << ")" << "+";
		}
		else if (poly[i+1] < 0) {
			cout << poly[i] << "*z^" << "(" << highest-i << ")" ;
		}
	}
	cout << endl;

}
	int degree()  {
	return deg;
}
	int highdeg() {
	return highest;
}
	void LaurentAdd(Laurent &A,Laurent &B){
	int ha=A.highdeg();
	int hb=B.highdeg();

	int lenA=A.degree()+1;
	int lenB=B.degree()+1;

	vector<T> coefA=A.poly;
	vector<T> coefB=B.poly;

	int lc,hc;

	if (ha > hb) {
		hc=ha;
		} else {
		hc=hb;
		}

	if (ha - A.deg < hb - B.deg ) {
		lc=ha-A.deg;
	} else {
		lc=ha-B.deg;
	}

	vector<T> coef_c;

	for (int i=0; i < hc-lc+1; i++) {
		coef_c.push_back(0);
	}

	for (int i=0; i < lenA; i++) {
		coef_c[hc-ha+i]=coef_c[hc-ha+i]+coefA[i];
	}

	for (int i=0; i < lenB; i++) {
		coef_c[hc-hb+i]=coef_c[hc-hb+i]+coefB[i];
	}



    setPoly(coef_c,hc);

}
	void LaurentSub(Laurent &A,Laurent &B){
	int ha=A.highdeg();
	int hb=B.highdeg();

	int lenA=A.degree()+1;
	int lenB=B.degree()+1;

	vector<T> coefA=A.poly;
	vector<T> coefB=B.poly;

	int lc,hc;

	if (ha > hb) {
		hc=ha;
		} else {
		hc=hb;
		}

	if (ha - A.deg < hb - B.deg ) {
		lc=ha-A.deg;
	} else {
		lc=ha-B.deg;
	}

	vector<T> coef_c;

	for (int i=0; i < hc-lc+1; i++) {
		coef_c.push_back(0);
	}

	for (int i=0; i < lenA; i++) {
		coef_c[hc-ha+i]=coef_c[hc-ha+i]+coefA[i];
	}

	for (int i=0; i < lenB; i++) {
		coef_c[hc-hb+i]=coef_c[hc-hb+i]-coefB[i];
	}



    setPoly(coef_c,hc);

}
   void LaurentMult(Laurent &A, Laurent &B) {
    int la=A.highdeg();
	int lb=B.highdeg();

	int lenA=A.degree()+1;
	int lenB=B.degree()+1;

	vector<T> coefA=A.poly;
	vector<T> coefB=B.poly;
	int degA,degB;
	vector<int> degAB;
	vector<T> coefAB;

	for (int i=0; i < lenA; i++) {
	    degA=i+la;
	    for (int j=0; j < lenB; j++) {
	        degB=j+lb;
	        degAB.push_back(degA+degB);
	        coefAB.push_back(coefA[i]*coefB[j]);

	    }

	}

	int min_deg= *min_element(degAB.begin(),degAB.end());
	int max_deg= *max_element(degAB.begin(),degAB.end());

	vector<T> coef_c;

	for (int i=0; i < abs(max_deg - min_deg)+1; i++) {
		coef_c.push_back(0);
	}

	for (int i=min_deg; i < max_deg+1;i++) {
	    for (int j=0; j < lenA*lenB;j++) {
	        if (i==degAB[j]) {
	            coef_c[i-min_deg]+=coefAB[j];
             }
        }
    }

	setPoly(coef_c,min_deg);


   }

   void zinv(Laurent &A) {
    int la=A.highdeg();

	//int lenA=A.degree()+1;

	vector<T> coefA=A.poly;

	la=(la-A.degree())*-1;

	reverse(coefA.begin(),coefA.end());

    setPoly(coefA,la);

   }

   bool isMono() {
    int la=highest;

	int lenA=deg+1;

	vector<T> coefA=poly;
	bool mono;

	if (lenA == 0) {
	    mono = false;
	}

	if (lenA == 1) {
	    if (coefA[0] != 0 || abs(coefA[0]) <= 1e-05) {
	        mono = true;

	    } else {
	        mono = false;

	    }

	}

	if (lenA > 1) {
	    int j=0;
	    for (int i=0; i < lenA; i++) {
	        if (coefA[i] == 0 || abs(coefA[i]) <= 1e-05) {
	            j++;
	        }

	    }

	    if (lenA==j+1) {
	        mono = true;

	    } else {
	        mono = false;
	    }

	}

    return mono;

}

void LaurentDiv(Laurent &A, Laurent &B, vector<Laurent<T> > &lcont) {
    int ha=A.highdeg();
	int hb=B.highdeg();

	int lenA=A.degree()+1;
	int lenB=B.degree()+1;

	int la=ha-lenA+1;
	int lb=hb-lenB+1;

	int hc = ha - hb;
	int lc = la -lb;

	int lenC = lenA - lenB;


	vector<T> coefA=A.poly;
	vector<T> coefB=B.poly;
	vector<T> coef_q1,coef_q2,coef_r1,coef_r2;

	if (lenC == 0) {
	    T temp1,temp2;
	    temp1=coefA[0]/coefB[0];
	    temp2=coefA[lenA-1]/coefB[lenB-1];

	    coef_q1.push_back(temp1);
	    coef_q2.push_back(temp2);

	    for (int i=0; i < lenA; i++) {
	        coef_r1.push_back(coefA[i]- temp1*coefB[i]);
	        coef_r2.push_back(coefA[i]- temp2*coefB[i]);

	    }
		Laurent<T> q1,q2,r1,r2;
		q1.setPoly(coef_q1,ha-hb);
		q2.setPoly(coef_q2,la-lb);
		r1.setPoly(coef_r1,ha);
		r2.setPoly(coef_r2,ha);

		//vector<Laurent<T> > lcont;
		lcont.push_back(q1);
		lcont.push_back(r1);
		lcont.push_back(q2);
		lcont.push_back(r2);



	}

    }

    virtual ~Laurent(){
    }
};

template <typename T>
class LaurentMat {
    Laurent<T> A,B,C,D;


public:
	LaurentMat() {

	}

	void setMat(Laurent<T> &AA, Laurent<T> &BB, Laurent<T> &CC, Laurent<T> &DD) {
	A=AA;
	B=BB;
	C=CC;
	D=DD;
	}

	void Det(Laurent<T> &oup) {
	    Laurent<T> tempMat1,tempMat2;
	    tempMat1.LaurentMult(A,D);
	    tempMat2.LaurentMult(C,B);
	    tempMat1.LaurentSub(tempMat1,tempMat2);
	    oup=tempMat1;

	}

    virtual ~LaurentMat() {

    }

};

template <typename T>
void Div(Laurent<T> &A, Laurent<T> &B, vector<Laurent<T> > &lcont) {
    int ha=A.highdeg();
	int hb=B.highdeg();

	int lenA=A.degree()+1;
	int lenB=B.degree()+1;

	int la=ha-lenA+1;
	int lb=hb-lenB+1;

	int hc = ha - hb;
	int lc = la -lb;

	int lenC = lenA - lenB;


	vector<T> coefA,coefB;
	A.getPoly(coefA);
	B.getPoly(coefB);
	
	if (lenC > 0) {
		vector<T> coef_q1,coef_q2,coef_r1,coef_r2;
		//vector<T> temp1,temp2;
		T t1,t2;
		coef_r1=coefA;
		
		for (int i=0; i < lenC+1; i++) {
			t1=coef_r1[i]/coefB[0];
			if (i == lenC) {
				t2=coef_r1[lenA-1]/coefB[lenB-1];
				coef_r2=coef_r1;
				int k=0;
			for (int j=lenA-lenB; j < lenA; j++) {
				coef_r2[j]=coef_r2[j]-t2*coefB[k];
				k++;
			}
			     coef_q2=coef_q1;
				 coef_q2.push_back(t2);
				
			}
			coef_q1.push_back(t1);
			int k=0;
			for (int j=i; j < i+lenB; j++) {
				coef_r1[j]=coef_r1[j]-t1*coefB[k];
				k++;
			}
			
		}
		    Laurent<T> q1,q2,r1,r2;
		    q1.setPoly(coef_q1,ha-hb);
		    q2.setPoly(coef_q2,ha-hb);
		    r1.setPoly(coef_r1,ha);
		    r2.setPoly(coef_r2,ha);

		//vector<Laurent<T> > lcont;
		    lcont.push_back(q1);
		    lcont.push_back(r1);
		    lcont.push_back(q2);
			lcont.push_back(r2);
			
			coef_q1.clear(),coef_q2.clear(),coef_r1.clear(),coef_r2.clear();
			coef_r1=coefA;
			
			for (int i=0; i < lenC+1; i++) {
			t1=coef_r1[lenA-i-1]/coefB[lenB-1];
			if (i == lenC) {
				t2=coef_r1[0]/coefB[0];
				coef_r2=coef_r1;
				int k=0;
			for (int j=0; j < lenB; j++) {
				coef_r2[j]=coef_r2[j]-t2*coefB[k];
				k++;
			}
			     coef_q2=coef_q1;
				 coef_q2.push_back(t2);
				
			}
			coef_q1.push_back(t1);
			int k=0;
			for (int j=lenC-i; j < lenA-i; j++) {
				coef_r1[j]=coef_r1[j]-t1*coefB[k];
				k++;
			}
			
		}
		    reverse(coef_q1.begin(),coef_q1.end());
			reverse(coef_q2.begin(),coef_q2.end());
			
			Laurent<T> q3,q4,r3,r4;
		    q3.setPoly(coef_q1,ha-hb);
		    q4.setPoly(coef_q2,(int)la-lb+coef_q1.size()-1);
		    r3.setPoly(coef_r1,ha);
		    r4.setPoly(coef_r2,ha);

		//vector<Laurent<T> > lcont;
		    lcont.push_back(q3);
		    lcont.push_back(r3);
		    lcont.push_back(q4);
			lcont.push_back(r4);
		
		
		
	} else if (lenC == 0) {
		vector<T> coef_q1,coef_q2,coef_r1,coef_r2;
	    T temp1,temp2;
	    temp1=coefA[0]/coefB[0];
	    temp2=coefA[lenA-1]/coefB[lenB-1];

	    coef_q1.push_back(temp1);
	    coef_q2.push_back(temp2);

	    for (int i=0; i < lenA; i++) {
	        coef_r1.push_back(coefA[i]- temp1*coefB[i]);
	        coef_r2.push_back(coefA[i]- temp2*coefB[i]);

	    }
		Laurent<T> q1,q2,r1,r2;
		q1.setPoly(coef_q1,ha-hb);
		q2.setPoly(coef_q2,la-lb);
		r1.setPoly(coef_r1,ha);
		r2.setPoly(coef_r2,ha);

		//vector<Laurent<T> > lcont;
		lcont.push_back(q1);
		lcont.push_back(r1);
		lcont.push_back(q2);
		lcont.push_back(r2);

	} else {
		vector<T> coef_q;
		coef_q.push_back(0);
		Laurent<T> q;
		q.setPoly(coef_q,0);

		//vector<Laurent<T> > lcont;
		lcont.push_back(q);
		lcont.push_back(B);
		
	}

    }
	

#endif /* LAURENT_H_ */


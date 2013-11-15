#ifndef ITERTEST_H
#define ITERTEST_H

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

template <typename T>
class fft_data{
public:	
	vector<T> re;
	vector<T> im;

};

template <typename T>
void inline twiddle(fft_data<T> &vec,int N,int radix){
	// Calculates twiddle factors for radix-2
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) PI2/N;
	vec.re.resize(N/radix,(T) 0.0);
	vec.im.resize(N/radix,(T) 0.0);
	vec.re[0] = (T) 1.0;
	
	for (int K = 1; K < N/radix; K++) {
		vec.re[K] = (T) cos(theta * K);
		vec.im[K] = (T) sin(theta * K);
	}

	
}

template <typename T>
void inline sh_radix5_dif(fft_data<T> &x,fft_data<T> &wl, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(5.0, (double)q);
	int Ls = L / 5;
	int r = n / L;
	
	T c1 = 0.30901699437;
	T c2 = -0.80901699437;
	T s1 = 0.95105651629;
	T s2 = 0.58778525229;
	
	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T br,bi,cr,ci,dr,di,er,ei;

	fft_data<T> y = x;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i;
	int lsr = Ls*r;
	
	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;
		
		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;
		
		wl4r = wl2r*wl2r - wl2i*wl2i;
		wl4i = 2.0*wl2r*wl2i;
		
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;

			tau0r = y.re[index1] + y.re[index4];
			tau0i = y.im[index1] + y.im[index4];
			
			tau1r = y.re[index2] + y.re[index3];
			tau1i = y.im[index2] + y.im[index3];
			
			tau2r = y.re[index1] - y.re[index4];
			tau2i = y.im[index1] - y.im[index4];
			
			tau3r = y.re[index2] - y.re[index3];
			tau3i = y.im[index2] - y.im[index3];
			
			tau4r = c1 * tau0r + c2 * tau1r;
			tau4i = c1 * tau0i + c2 * tau1i;
			
			tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
			tau5i = sgn * ( s1 * tau2i + s2 * tau3i);
			
			br = y.re[index] + tau4r + tau5i;
			bi = y.im[index] + tau4i - tau5r;
			
			er = y.re[index] + tau4r - tau5i;
			ei = y.im[index] + tau4i + tau5r;
			
			tau4r = c2 * tau0r + c1 * tau1r;
			tau4i = c2 * tau0i + c1 * tau1i;
			
			tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
			tau5i = sgn * ( s2 * tau2i - s1 * tau3i);
			
			cr = y.re[index] + tau4r + tau5i;
			ci = y.im[index] + tau4i - tau5r;
			
			dr = y.re[index] + tau4r - tau5i;
			di = y.im[index] + tau4i + tau5r;

			int indexo = k*Ls+j;
			int indexo1 = indexo+lsr;
			int indexo2 = indexo1+lsr;
			int indexo3 = indexo2+lsr;
			int indexo4 = indexo3+lsr;

			x.re[indexo]= y.re[index] + tau0r + tau1r;
			x.im[indexo]= y.im[index] + tau0i + tau1i;

			x.re[indexo1] = wlr*br - wli*bi;
			x.im[indexo1] = wlr*bi + wli*br;
			
			x.re[indexo2] = wl2r*cr - wl2i*ci;
			x.im[indexo2] = wl2r*ci + wl2i*cr;

			x.re[indexo3] = wl3r*dr - wl3i*di;
			x.im[indexo3] = wl3r*di + wl3i*dr;

			x.re[indexo4] = wl4r*er - wl4i*ei;
			x.im[indexo4] = wl4r*ei + wl4i*er;
			
			
		}

	}
	
}

template <typename T>
void inline fftsh_radix5_dif(fft_data<T> &data,int sgn, unsigned int N) {
	//unsigned int len = data.re.size(); 
	
	int num = (int) ceil(log10(static_cast<double>(N))/log10(5.0));

	//indrev(data,index);
	fft_data<T> twi;
	
	twiddle(twi,N,5);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	for (int i=num; i > 0; i--) {
    sh_radix5_dif(data,twi,i,sgn);
	
	
    }
	
	
}


#endif
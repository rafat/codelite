
#include "itertest.h"

using namespace std;



int main(int argc, char **argv)
{
	int N = 25;
	//vector<complex<double> > sig1;
	fft_data<double> sig1;
	for (int i =0; i < N; i++){
	//sig1.push_back(complex<double>((double)1.0, 0.0));
		//sig2.re.push_back((double) i);
		//sig2.im.push_back((double) i+2);
		sig1.re.push_back((double) 1);
		sig1.im.push_back((double) 0);
		//cout << real(sig1[i]) << " " << imag(sig1[i]) << endl;
	}

	fftsh_radix5_dif(sig1,1,N);

	for (int i =0; i < N; i++){
		cout << sig1.re[i] << " " << sig1.im[i] << endl;
	}
     
	return 0;

}
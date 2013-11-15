
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

int main() {
  int n = 4;      
  static const double pi=3.14159265358979323846;
  size_t n2=8;
  size_t m, coeff;
  double angle, xn2;
  vector<double> bk(100), bkf(100), work(100);
  double pibyn=pi/n;
  

/* initialize b_k */
  bk[0] = 1;
  bk[1] = 0;

  coeff=0;
  for (m=1; m<n; ++m)
    {
    coeff+=2*m-1;
    if (coeff>=2*n) coeff-=2*n;
    angle = pibyn*coeff;
    bk[2*m] = cos(angle);
    bk[2*m+1] = sin(angle);
	
    }

/* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
  xn2 = 1.0;
  bkf[0] = bk[0]*xn2;
  bkf[1] = bk[1]*xn2;
  for (m=2; m<2*n; m+=2)
    {
    bkf[m]   = bkf[2*n2-m]   = bk[m]   *xn2;
    bkf[m+1] = bkf[2*n2-m+1] = bk[m+1] *xn2;
	
    }
  for (m=2*n;m<=(2*n2-2*n+1);++m)
    bkf[m]=0.;

  for (int i = 0; i < 100; i++)
	  cout << bkf[2*i] << " " << bkf[2*i+1] << endl;
            
            //gnudwtplot(J);
            return 0;
}
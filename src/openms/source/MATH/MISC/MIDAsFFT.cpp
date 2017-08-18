#include <OpenMS/MATH/MISC/MIDAsFFT.h>

namespace OpenMS
{

  void fft(double *Data, int nn)
  {
    int isign = -1;
    unsigned long i, j, m, n, mmax, istep;
    double wr, wpr, wpi, wi, theta;
    double wtemp, tempr, tempi;
    double one_pi = acos(-1);
    double two_pi= 2*one_pi;

    /* Perform bit reversal of Data[] */
    n = nn << 1;
    j=1;
    for (i=1; i<n; i+=2)
    {
      if (j > i)
      {
        wtemp = Data[i];
        Data[i] = Data[j];
        Data[j] = wtemp;
        wtemp = Data[i+1];
        Data[i+1] = Data[j+1];
        Data[j+1] = wtemp;
      }
      m = n >> 1;
      while (m >= 2 && j > m)
      {
        j -= m;
        m >>= 1;
      }
      j += m;
    }
 
    /* Perform Danielson-Lanczos section of FFT */
    n = nn << 1;
    mmax = 2;
    while (n > mmax)  /* Loop executed log(2)nn times */
    {
      istep = mmax << 1;
      theta = isign * (two_pi/mmax);  
      wtemp = sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for (m=1; m<mmax; m+=2)
      {
        for (i=m; i<=n; i+=istep)
        {
          j = i+mmax;                      

          tempr = wr*Data[j]-wi*Data[j+1];
          tempi = wr*Data[j+1]+wi*Data[j];
          Data[j] = Data[i]-tempr;
          Data[j+1] = Data[i+1]-tempi;
          Data[i] += tempr;
          Data[i+1] += tempi;
        }
        wr = (wtemp=wr)*wpr-wi*wpi+wr;
        wi = wi*wpr+wtemp*wpi+wi;
      }
      mmax = istep;
    }

  }

}

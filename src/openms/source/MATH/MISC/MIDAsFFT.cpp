#include <OpenMS/MATH/MISC/MIDAsFFT.h>


#define INVERSE
using namespace std;

namespace OpenMS
{

  void fft(double *Data, int nn)
  {
    unsigned long i, j, m, n, mmax;

    /* Perform bit reversal of Data[] */
    n = nn << 1;
    j=1;
    for (i=1; i<n; i+=2)
    {
      if (j > i)
      {
        swap(Data[i], Data[j]);
        swap(Data[i+1], Data[j+1]);
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
      unsigned long istep = mmax << 1;
#ifdef INVERSE
      double theta = -((2 * Constants::PI) / mmax);
#else
      double theta = ((2 * Constants::PI) / mmax);
#endif
      double wtemp = sin(0.5 * theta);
      double wpr = -2.0 * wtemp * wtemp;
      double wpi = sin(theta);
      double wr = 1.0;
      double wi = 0.0;
      for (m=1; m<mmax; m+=2)
      {
        for (i=m; i<=n; i+=istep)
        {
          double tempr, tempi;
          j = i+mmax;                      
          tempr = wr * Data[j] - wi * Data[j + 1];
          tempi = wr * Data[j + 1] + wi * Data[j];
          Data[j] = Data[i] - tempr;
          Data[j + 1] = Data[i + 1] - tempi;
          Data[i] += tempr;
          Data[i + 1] += tempi;
        }
        wtemp = wr;
        wr += wtemp * wpr - wi * wpi;
        wi += (wi * wpr) + (wtemp * wpi);
      }
      mmax = istep;
    }

  }

}

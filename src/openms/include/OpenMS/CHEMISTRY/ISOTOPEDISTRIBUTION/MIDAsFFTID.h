 

#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_MIDASFFTID_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_MIDASFFTID_H



// defines required for kissfft
#define kiss_fft_scalar double
#define INVERSE true

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>
#include <kiss_fft.h>


namespace OpenMS
{
  class OPENMS_DLLAPI MIDAsFFTID : public MIDAs
  {
 public:
    typedef kiss_fft_cpx fft_complex;
    typedef std::vector<fft_complex> FFT_Spectrum;
    typedef struct {double mean; double variance;} Stats;
    MIDAsFFTID(EmpiricalFormula&, double);
    void init();
    void run();
 private:
    FFT_Spectrum input_, output_;
   
    double cutoff_amplitude_factor_;
    double average_mass_;

    double delta_;
    double mass_range_;
    //void (double);
    Stats formulaMeanAndVariance(double resolution = 1.0);
   
  };
}
#endif

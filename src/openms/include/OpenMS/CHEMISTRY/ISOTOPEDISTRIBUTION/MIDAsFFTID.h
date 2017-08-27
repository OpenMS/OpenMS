 

#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_MIDASFFTID_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_MIDASFFTID_H



// defines required for kissfft
#define kiss_fft_scalar double

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
    MIDAsFFTID(double);
    void init(const EmpiricalFormula&);
    void run(const EmpiricalFormula&);
 private:
    FFT_Spectrum input_, output_;
    double cutoff_amplitude_factor_;
    double average_mass_;
    double delta_;
    double mass_range_;
    Stats formulaMeanAndVariance(const EmpiricalFormula&, double resolution = 1.0);
   
  };
}

#endif

#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_ECIPEX_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_ECIPEX_H

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Base.h>


namespace OpenMS
{
  /**
        @ingroup Chemistry

        @brief Isotope distribution class

        Holds an isotope distribution with the weight value and according
        probability. Distribution can be add using the '+' or '+=' operators.

        The most important value which should be set is the max isotope value.
        This value can be set using the setMaxIsotope method. It is an upper
        bound for the number of isotopes which are calculated. E.g. if it is set
        to 3, only the first three isotopes, Monoisotopic mass, +1 and +2 are
        calculated.
        By default all possible isotopes are calculated, which leads to a large
        number of values, if the mass value is large!
    */

  class OPENMS_DLLAPI Ecipex : public IsotopeDistribution
  {
  public:
    typedef ContainerType Spectrum;
    using IsotopeDistribution::operator[];
    Ecipex(EmpiricalFormula& formula, double threshold, double fft_threshold);
    Ecipex();
    Ecipex(const IsotopeDistribution& isotope_distribution);
    
    void sortAndNormalize();
    void computeIsotopePattern(double threshold, double fft_threshold);
    void run();
    Ecipex elementIsotopePattern(const Spectrum& iso_pattern, UInt size, double fft_threshold);
    void convolve(IsotopeDistribution& spectrum, double threshold);
  private:
    EmpiricalFormula formula_;
    double fft_threshold_;
    double threshold_;
  };


}




#endif // OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_ECIPEX_H

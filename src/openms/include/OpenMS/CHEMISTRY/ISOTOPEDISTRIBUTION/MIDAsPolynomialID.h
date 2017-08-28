
#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_MIDASPOLYNOMIALID_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_MIDASPOLYNOMIALID_H

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>
  
namespace OpenMS
{

  class OPENMS_DLLAPI MIDAsPolynomialID : public MIDAs
  {
 public:

    MIDAsPolynomialID(double resolution, double probability_cutoff);
    void run(const EmpiricalFormula&);
    
 private:
    Polynomial generatePolynomial(const Element&, const SignedSize);
    void multiplyPolynomials(Polynomial&, Polynomial&);
    void merge_polynomial(Polynomial&);
    
  };
}

#endif

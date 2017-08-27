
#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_MIDASPOLYNOMIALID_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_MIDASPOLYNOMIALID_H

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>
  
namespace OpenMS
{

  class OPENMS_DLLAPI MIDAsPolynomialID : public MIDAs
  {
 public:

    MIDAsPolynomialID(double);
    void run(const EmpiricalFormula&);
    
 private:
    Polynomial generatePolynomial(const Element&, const SignedSize);
    void multiplyPolynomials(Polynomial&, Polynomial&);
    void merge_polynomial(Polynomial&);
    void dumpID(Polynomial&);
    
    //Polynomial fgid;
    double fine_resolution;
    
    double mw_resolution;
    double resolution;
    double min_resolution;
    
  };
}

#endif

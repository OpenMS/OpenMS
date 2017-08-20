
#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_MIDASPOLYNOMIALID_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_MIDASPOLYNOMIALID_H

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternSolver.h>
  
namespace OpenMS
{

  class OPENMS_DLLAPI MIDAsPolynomialID : public MIDAs
  {
 public:

    MIDAsPolynomialID(EmpiricalFormula&, double);
    void run();
    
 private:
    Polynomial generatePolynomial(const Element&, const SignedSize);
    double lightest_mass();
    void multiplyPolynomials(Polynomial&, Polynomial&);
    void merge_polynomial(Polynomial&);
    void dumpID(Polynomial&);
    double fact_ln(UInt);
    
    //Polynomial fgid;
    double fine_resolution;

    double lighter_isotope;
    
    double mw_resolution;
    double resolution;
    double min_resolution;
    
  };
}

#endif

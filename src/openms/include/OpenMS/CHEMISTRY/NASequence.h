#ifndef OPENMS_CHEMISTRY_NASEQUENCE_H
#define OPENMS_CHEMISTRY_NASEQUENCE_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <vector>
#include <iosfwd>

namespace OpenMS
{
    class OPENMS_DLLAPI NASequence
  {
    public:
       NASequence(); //default constructor
       NASequence(const String & rhs); // copy constructor
       NASequence& operator=(const NASequence& rhs); //assignment operator
       bool operator==(const NASequence& rhs) const;
       virtual ~NASequence(); //destructor
       void setSequence(const String & s);
       String getSequence() const;
       size_t size() const;
       double getMonoWeight(Residue::ResidueType type, Int charge) const;
       NASequence getPrefix(Size index) const;
       NASequence getSuffix(Size index) const;
       bool empty() const;
       EmpiricalFormula getFormula(Residue::ResidueType type, Int charge) const;
    private:
        String s_;
    };
OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const NASequence& nucleotide);

OPENMS_DLLAPI std::istream& operator>>(std::istream& os, const NASequence& nucleotide);
}
#endif // OPENMS_CHEMISTRY_NASEQUENCE_H

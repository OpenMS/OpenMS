// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Timo Sachsenberg $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_ELEMENTDB_H
#define OPENMS_CHEMISTRY_ELEMENTDB_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Container.h>
#include <OpenMS/CHEMISTRY/Element.h>

namespace OpenMS
{

  /** @ingroup Chemistry

          @brief Stores elements

      The elements weights (in the default file) are taken from
      "Isotopic Compositions of the Elements 1997", Pure Appl. Chem., 70(1), 217-235, 1998.
      (http://www.iupac.org/reports/1998/7001rosman/)

      The isotope distributions (in the default file) are taken from
          "Atomic weights of the elements. Review 2000" (IUPAC Technical Report)
          Pure Appl. Chem., 2003, Vol. 75, No. 6, pp. 683-799
          doi:10.1351/pac200375060683

          This singleton stores all elements. The elements are taken from the publications given
          above and are stored in share/OpenMS/CHEMISTRY/Elements.xml.

          Specific isotopes of elements can be accessed by writing the atomic number of the isotope
          in brackets followed by the element name, e.g. "(2)H" for deuterium.

    @improvement include exact mass values for the isotopes (done) and update IsotopeDistribution (Andreas)
          @improvement add exact isotope distribution based on exact isotope values (Andreas)
*/

  class OPENMS_DLLAPI ElementDB
  {
public:

    /** @name Accessors
    */
    //@{
    /// returns a pointer to the singleton instance of the element db
    inline static const ElementDB * getInstance()
    {
      static ElementDB * db_ = nullptr;
      if (db_ == nullptr)
      {
        db_ = new ElementDB;
      }
      return db_;
    }

    /// returns a hashmap that contains names mapped to pointers to the elements
    const Map<String, const Element *> & getNames() const;

    /// returns a hashmap that contains symbols mapped to pointers to the elements
    const Map<String, const Element *> & getSymbols() const;

    /// returns a hashmap that contains atomic numbers mapped to pointers of the elements
    const Map<UInt, const Element *> & getAtomicNumbers() const;

    /** returns a pointer to the element with name or symbol given in parameter name;
        *	if no element exists with that name or symbol 0 is returned
        *	@param name: name or symbol of the element
    */
    const Element * getElement(const String & name) const;

    /// returns a pointer to the element of atomic number; if no element is found 0 is returned
    const Element * getElement(UInt atomic_number) const;

    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if the db contains an element with the given name
    bool hasElement(const String & name) const;

    /// returns true if the db contains an element with the given atomic_number
    bool hasElement(UInt atomic_number) const;
    //@}

protected:

    /*_ parses a Histogram given as a OpenMS String and return the distribution

            @throw throws exception ParseError
     */
    IsotopeDistribution parseIsotopeDistribution_(const Map<UInt, double>& Z_to_abundance, const Map<UInt, double>& Z_to_mass);

    /*_ calculates the average weight based on isotope abundance and mass
     */
    double calculateAvgWeight_(const Map<UInt, double> & Z_to_abundance, const Map<UInt, double> & Z_to_mass);

    /*_ calculates the mono weight based on the smallest isotope mass
     */
    double calculateMonoWeight_(const Map<UInt, double> & Z_to_mass);

    /*_ read elements from a XML file, formatted as a Param file.

            @throw throws ParseError if the file cannot be parsed
            @throw throws FileNotFound if the file could not be found
     */
    void readFromFile_(const String & file_name);

    /*_ resets all containers
     */
    void clear_();

    Map<String, const Element *> names_;

    Map<String, const Element *> symbols_;

    Map<UInt, const Element *> atomic_numbers_;

private:

    ElementDB();

    ElementDB(const ElementDB & db);

    ElementDB & operator=(const ElementDB & db);

    virtual ~ElementDB();
  };

} // namespace OpenMS
#endif

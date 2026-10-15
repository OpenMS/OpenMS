// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Timo Sachsenberg, Chris Bielow, Jang Jang Jin$
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

#include <map>
#include <memory>
#include <unordered_map>
#include <string>

namespace OpenMS
{
  class Element;

  /** @ingroup Chemistry

      @brief Singleton that stores elements and isotopes.

      The elements weights (in the default file) are taken from
      "Isotopic Compositions of the Elements 1997", Pure Appl. Chem., 70(1), 217-235, 1998.
      (http://www.iupac.org/reports/1998/7001rosman/)

      The isotope distributions (in the default file) are taken from
          "Atomic weights of the elements. Review 2000" (IUPAC Technical Report)
          Pure Appl. Chem., 2003, Vol. 75, No. 6, pp. 683-799
          doi:10.1351/pac200375060683

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
    /// returns a pointer to the (immutable) singleton instance of the element db
    /// This is thread safe upon first and subsequent calls.
    static const ElementDB* getInstance();

    /// returns a hashmap that contains names mapped to pointers to the elements
    const std::unordered_map<std::string, const Element*>& getNames() const;

    /// returns a hashmap that contains symbols mapped to pointers to the elements
    const std::unordered_map<std::string, const Element*>& getSymbols() const;

    /// returns a hashmap that contains atomic numbers mapped to pointers of the elements
    const std::unordered_map<unsigned int, const Element*>& getAtomicNumbers() const;

    /** returns a pointer to the element with name or symbol given in parameter name;
        *	if no element exists with that name or symbol 0 is returned
        *	@param[in] name name or symbol of the element
    */
    const Element* getElement(const std::string& name) const;

    /// returns a pointer to the element of atomic number; if no element is found 0 is returned
    const Element* getElement(unsigned int atomic_number) const;
    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if the db contains an element with the given name
    bool hasElement(const std::string& name) const;

    /// returns true if the db contains an element with the given atomic_number
    bool hasElement(unsigned int atomic_number) const;
    //@}

protected:

    /** parses a isotope distribution of abundances and masses

    **/
    IsotopeDistribution parseIsotopeDistribution_(const std::map<unsigned int, double>& abundance, const std::map<unsigned int, double>& mass);

    /** calculates the average weight based on isotope abundance and mass
     **/
    double calculateAvgWeight_(const std::map<unsigned int, double>& abundance, const std::map<unsigned int, double>& mass);

    /**_ calculates the mono weight based on the most abundant isotope 
     **/
    double calculateMonoWeight_(const std::map<unsigned int, double>& abundance, const std::map<unsigned int, double>& mass);

	  /// constructs element objects
    void storeElements_();

    /// build element objects from given abundances, masses, name, symbol, and atomic number
    void buildElement_(const std::string& name, const std::string& symbol, const unsigned int an, const std::map<unsigned int, double>& abundance, const std::map<unsigned int, double>& mass);

    /// add element objects to documentation maps
    void addElementToMaps_(const std::string& name, const std::string& symbol, const unsigned int an, std::unique_ptr<const Element> e);

    /// constructs isotope objects
    void storeIsotopes_(const std::string& name, const std::string& symbol, const unsigned int an, const std::map<unsigned int, double>& Z_to_mass, const IsotopeDistribution& isotopes);

    /**_ resets all containers
    **/
    void clear_();

    std::unordered_map<std::string, const Element*> names_;

    std::unordered_map<std::string, const Element*> symbols_;

    std::unordered_map<unsigned int, const Element*> atomic_numbers_;

private:
    ElementDB();
    ~ElementDB();
    ElementDB(const ElementDB& db) = delete;
    ElementDB(const ElementDB&& db) = delete;
    ElementDB& operator=(const ElementDB& db) = delete;

  };

} // namespace OpenMS

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Timo Sachsenberg, Chris Bielow, Jang Jang Jin$
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
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
    /// returns a pointer to the singleton instance of the element db
    /// This is thread safe upon first and subsequent calls.
    static ElementDB* getInstance();

    /// returns a hashmap that contains names mapped to pointers to the elements
    const std::unordered_map<std::string, const Element*>& getNames() const;

    /// returns a hashmap that contains symbols mapped to pointers to the elements
    const std::unordered_map<std::string, const Element*>& getSymbols() const;

    /// returns a hashmap that contains atomic numbers mapped to pointers of the elements
    const std::unordered_map<unsigned int, const Element*>& getAtomicNumbers() const;

    /** returns a pointer to the element with name or symbol given in parameter name;
        *	if no element exists with that name or symbol 0 is returned
        *	@param name: name or symbol of the element
    */
    const Element* getElement(const std::string& name) const;

    /// returns a pointer to the element of atomic number; if no element is found 0 is returned
    const Element* getElement(unsigned int atomic_number) const;

    /** Adds or replaces a new element to the database
     *
     * Adds a new element (or replaces an existing one if @em replace_existing is true). 
     *
     * @param name Common name of the element
     * @param symbol Element symbol (one or two letter)
     * @param an Atomic number (number of protons)
     * @param abundance List of abundances for each isotope (e.g. {{12u, 0.9893}, {13u, 0.0107}} for Carbon)
     * @param mass List of masses for each isotope (e.g. {{12u, 12.0}, {13u, 13.003355}} for Carbon)
     * @param replace_existing If the element must be replaced (i.e. is not new), either allow that (=true), or throw an exception (=false).
     *
     * @throw Exception::IllegalArgument if element already exists in DB, but @p replace_existing is false
     *
     * @note Do not use this function inside parallel code as it modifies a singleton that is shared between threads.
    */
    void addElement(const std::string& name,
                    const std::string& symbol,
                    const unsigned int an,
                    const std::map<unsigned int, double>& abundance,
                    const std::map<unsigned int, double>& mass,
                    bool replace_existing);
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

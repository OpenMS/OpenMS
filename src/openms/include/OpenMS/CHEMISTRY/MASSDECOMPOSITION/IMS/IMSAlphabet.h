// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#pragma once

#include <vector>
#include <string>
#include <iosfwd>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetParser.h>

namespace OpenMS
{

  namespace ims
  {

    /**
      @brief Holds an indexed list of bio-chemical elements used by mass-decomposition algorithms.

      Stores @c IMSElement instances and exposes them in two complementary ways:
      lookup by index is constant-time, lookup by element name is linear-time.
      Because @c IMSAlphabet keeps the full element representation, prefer the
      lighter @c ims::Weights when only the mass values are needed in
      performance-critical code.

      An alphabet can be sorted either by element name or by mass; the file
      loaders always leave the alphabet sorted by mass, so files do not need
      to be pre-sorted.

      To populate an alphabet from a flat text file with the default format,
      call @c load(const std::string&). To support a different format, derive
      a parser from @c IMSAlphabetParser and pass it to the two-argument
      @c load() overload. Both overloads report a missing or unreadable file
      by throwing @c Exception::IOException.

      @ingroup Analysis_DeNovo
    */
    class OPENMS_DLLAPI IMSAlphabet
    {

public:
      typedef IMSElement element_type;
      typedef element_type::mass_type mass_type;
      typedef element_type::name_type name_type;
      typedef std::vector<element_type> container;
      typedef container::size_type size_type;
      typedef container::iterator iterator;
      typedef container::const_iterator const_iterator;
      typedef std::vector<name_type> name_container;
      typedef name_container::iterator name_iterator;
      typedef name_container::const_iterator const_name_iterator;
      typedef std::vector<mass_type> mass_container;
      typedef mass_container::iterator mass_iterator;
      typedef mass_container::const_iterator const_mass_iterator;
      typedef std::vector<mass_type> masses_type;

      /**
        Empty constructor.
      */
      IMSAlphabet() {}


      /**
        Constructor with elements.

        @param[in] elements Elements to be set
      */
      explicit IMSAlphabet(const container & elements) :
        elements_(elements)
      {}


      /**
        Copy constructor.

        @param[in] alphabet Alphabet whose elements are copied.
      */
      IMSAlphabet(const IMSAlphabet & alphabet) :
        elements_(alphabet.elements_)
      {}

      /**
        Returns the alphabet size.

        @return The size of alphabet.
      */
      size_type size() const
      {
        return elements_.size();
      }

      /**
        Gets the element with index @c index.
        @note Operation takes constant time.

        @param[in] index of the element
        @return Element with the given index in alphabet
      */
      const element_type & getElement(size_type index) const
      {
        return elements_[index];
      }

      /**
        Overwrites an element in the alphabet with the @c name with a new element constructed
        from the given name @c name and mass @c mass.
        If the parameter @c forced is set to true, a new element will be appended to the alphabet
        in the case the alphabet contains no element with the name @c name.

        @param[in] name The name of the element that should be replaced in (or appended to) the alphabet.
        @param[in] mass The new mass of the element in the alphabet.
        @param[in] forced Indicates whether a new element should be created (if set to @c true) if there is no element with the name @c name or not (if set to @c false).
      */
      void setElement(const name_type & name, mass_type mass, bool forced = false);

      /**
        Removes the element with name @c name from the alphabet.

        @param[in] name The name of the element to be removed from the alphabet.
        @return A boolean indicating whether an element was removed (@c true) or not (@c false).
      */
      bool erase(const name_type & name);

      /**
        Gets the element with the symbol @name. If there is
        no such element, throws @c Exception::InvalidValue.

        @param[in] name Name of the element.
        @return Element with the given name, or if there are no such element
        @throws Exception::InvalidValue.
      */
      const element_type & getElement(const name_type & name) const;

      /**
        Gets the symbol of the element with an index @c index in alphabet.

        @param[in] index of the element.
        @return Name of the element.
      */
      const name_type & getName(size_type index) const;

      /**
        Gets mono isotopic mass of the element with the symbol @c name.
        If there is no such element, throws an @c Exception::InvalidValue.

        @param[in] name Symbol of the element.
        @return Mass of the element, or if there are no element
        @throws Exception::InvalidValue.
        @see getMass(size_type index)
      */
      mass_type getMass(const name_type & name) const;

      /**
        Gets mass of the element with an index @c index in alphabet.

        @param[in] index Index of the element.
        @return Mass of the element.
        @see getMass(const std::string& name)
      */
      mass_type getMass(size_type index) const;

      /**
        Gets masses of elements isotopes given by @c isotope_index.

        @param[in] isotope_index Index of isotope
        @return Masses of elements isotopes with the given index.
      */
      masses_type getMasses(size_type isotope_index = 0) const;

      /**
        Gets average masses of elements.

        @return Average masses of elements.
      */
      masses_type getAverageMasses() const;

      /**
        Returns true if there is an element with symbol
        @c name in the alphabet, false - otherwise.

        @return True, if there is an element with symbol
                @c name, false - otherwise.
      */
      bool hasName(const name_type & name) const;

      /**
        Adds a new element with name @c name and mass @c value
        to the alphabet.

        @param[in] name Name of the element to be added.
        @param[in] value Mass of the element to be added.

        @see push_back(const element_type&)
      */
      void push_back(const name_type & name, mass_type value)
      {
        push_back(element_type(name, value));
      }

      /**
        Adds a new element @c element to the alphabet.

        @param[in] element The @c Element to be added.
      */
      void push_back(const element_type & element)
      {
        elements_.push_back(element);
      }

      /**
        Clears the alphabet data.
      */
      void clear()
      {
        elements_.clear();
      }

      /**
        Sorts the alphabet by names.

        @see sortByValues()
      */
      virtual void sortByNames();


      /**
        Sorts the alphabet by mass values.

        @see sortByNames()
      */
      virtual void sortByValues();


      /**
        Replaces the alphabet's contents with the elements read from the
        text file @p fname (default OpenMS alphabet text format), sorted
        by mass.

        @param[in] fname File to read.
        @throws Exception::IOException if @p fname cannot be opened.

        @see load(const std::string&, IMSAlphabetParser<>&)
      */
      virtual void load(const std::string & fname);


      /**
        Replaces the alphabet's contents with the elements read from
        @p fname through @p parser, sorted by mass.

        @param[in]     fname  File to read.
        @param[in,out] parser Parser used to read the file.
        @throws Exception::IOException if @p fname cannot be opened by @p parser.

        @see load(const std::string&)
        @see IMSAlphabetParser
      */
      virtual void load(const std::string & fname, IMSAlphabetParser<> & parser);


      /**
        Default destructor.
      */
      virtual ~IMSAlphabet() {}

private:
      /**
        Elements of the alphabet.
      */
      container elements_;

      /**
        @brief Private class-functor to sort out elements in mass ascending order.
      */
      class OPENMS_DLLAPI MassSortingCriteria_
      {
public:
        bool operator()(const element_type & el1,
                        const element_type & el2) const
        {
          return el1.getMass() < el2.getMass();
        }

      };

    };

    /**
      Prints alphabet to the stream @c os.

      @param[out] os Output stream to which alphabet is written
      @param[in] alphabet Alphabet to be written.
    */
    OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const IMSAlphabet & alphabet);

  } // namespace ims

} // namespace OpenMS


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
      @brief Holds an indexed list of bio-chemical elements.

      Presents an indexed list of bio-chemical elements of type (or derived from
      type) @c Element. Due to indexed structure @c Alphabet can be used similar
      to @c std::vector, for example to add a new element to @c Alphabet function
      @c push_back(element_type) can be used. Elements or their properties (such
      as element's mass) can be accessed by index in a constant time. On the other
      hand accessing elements by their names takes linear time. Due to this and
      also the fact that @c Alphabet is 'heavy-weighted' (consisting of
      @c Element -s or their derivatives where the depth of derivation as well is
      undefined resulting in possibly 'heavy' access operations) it is recommended
      not use @c Alphabet directly in operations where fast access to
      @c Element 's properties is required. Instead consider to use
      'light-weighted' equivalents, such as @c Weights.

      Elements in @c Alphabet can be sorted by the @c Element 's properties:
      sequence and mass. When alphabet's data is loaded from file it is
      automatically sorted by mass. To load data from file default function
      @c load(str::string& fname) can be used. Then elements have to be stored
      in a flat file @c fname in a predefined format. To read more on this format,
      please, @see AlphabetParser. If one wants to load data stored differently or
      in its own file format (i.e. xml) one has to define a new parser derived from
      @c AlphabetParser and pass its pointer together with the file name to function
      @c load(const std::string& fname, AlphabetParser<>* parser). If there is any
      error happened while loading data, @c IOException will be thrown.
    *
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

        @param elements Elements to be set
      */
      explicit IMSAlphabet(const container & elements) :
        elements_(elements)
      {}


      /**
        Copy constructor.

        @param alphabet Alphabet to be assigned
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

        @param index of the element
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

        @param name The name of the element that should be replaced in (or appended to) the alphabet.
        @param mass The new mass of the element in the alphabet.
        @param forced Indicates whether a new element should be created (if set to @c true) if there is no element with the name @c name or not (if set to @c false).
      */
      void setElement(const name_type & name, mass_type mass, bool forced = false);

      /**
        Removes the element with name @c name from the alphabet.

        @param name The name of the element to be removed from the alphabet.
        @return A boolean indicating whether an element was removed (@c true) or not (@c false).
      */
      bool erase(const name_type & name);

      /**
        Gets the element with the symbol @name. If there is
        no such element, throws @c Exception::InvalidValue.

        @param name Name of the element.
        @return Element with the given name, or if there are no such element
        @throws Exception::InvalidValue.
      */
      const element_type & getElement(const name_type & name) const;

      /**
        Gets the symbol of the element with an index @c index in alphabet.

        @param index of the element.
        @return Name of the element.
      */
      const name_type & getName(size_type index) const;

      /**
        Gets mono isotopic mass of the element with the symbol @c name.
        If there is no such element, throws an @c Exception::InvalidValue.

        @param name Symbol of the element.
        @return Mass of the element, or if there are no element
        @throws Exception::InvalidValue.
        @see getMass(size_type index)
      */
      mass_type getMass(const name_type & name) const;

      /**
        Gets mass of the element with an index @c index in alphabet.

        @param index Index of the element.
        @return Mass of the element.
        @see getMass(const std::string& name)
      */
      mass_type getMass(size_type index) const;

      /**
        Gets masses of elements isotopes given by @c isotope_index.

        @param isotope_index Index of isotope
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

        @param name Name of the element to be added.
        @param value Mass of the element to be added.

        @see push_back(const element_type&)
      */
      void push_back(const name_type & name, mass_type value)
      {
        push_back(element_type(name, value));
      }

      /**
        Adds a new element @c element to the alphabet.

        @param element The @c Element to be added.
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
        Loads the alphabet data from the file @c fname using the default
        parser. If there is no file @c fname, throws an @c IOException.

        @param fname The file name to be loaded.
        @throws Exception::IOException

        @see load(const std::string& fname, AlphabetParser<>* parser)
      */
      virtual void load(const std::string & fname);


      /**
        Loads the alphabet data from the file @c fname using @c parser.
        If there is no file @c fname found, throws an @c IOException.

        @param fname File name to be loaded.
        @param parser Parser to be used by loading.
        @throws Exception::IOException

        @see load(const std::string& fname)
        @see AlphabetParser
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

      @param os Output stream to which alphabet is written
      @param alphabet Alphabet to be written.
    */
    OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const IMSAlphabet & alphabet);

  } // namespace ims

} // namespace OpenMS


// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#include <functional>
#include <algorithm>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetTextParser.h>

namespace OpenMS::ims
{

    const IMSAlphabet::name_type & IMSAlphabet::getName(size_type index) const
    {
      return getElement(index).getName();
    }

    IMSAlphabet::mass_type IMSAlphabet::getMass(size_type index) const
    {
      return getElement(index).getMass();
    }

    IMSAlphabet::mass_type IMSAlphabet::getMass(const name_type & name) const
    {
      return getElement(name).getMass();
    }

    bool IMSAlphabet::hasName(const name_type & name) const
    {
      return std::find_if(elements_.begin(), elements_.end(),
                          [&name](const element_type& e) { return e.getName() == name; })
             != elements_.end();
    }

    const IMSAlphabet::element_type & IMSAlphabet::getElement(const name_type & name) const
    {
      for (const IMSElement& cit : elements_)
      {
        if (cit.getName() == name)
        {
          return cit;
        }
      }
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name + " was not found in IMSAlphabet!", String(name));
    }

    void IMSAlphabet::setElement(const name_type & name, mass_type mass, bool forced)
    {
      bool found = false;
      for (size_type i = 0; i < elements_.size(); ++i)
      {
        if (name == elements_[i].getName())
        {
          element_type element(name, mass);
          elements_[i] = element;
          found = true;
          break;
        }
      }
      if (!found && forced)
      {
        this->push_back(name, mass);
      }
    }

    bool IMSAlphabet::erase(const name_type & name)
    {
      bool found = false;
      for (iterator it = elements_.begin(); it != elements_.end(); ++it)
      {
        if (it->getName() == name)
        {
          elements_.erase(it);
          found = true;
          break;
        }
      }
      return found;
    }

    IMSAlphabet::masses_type IMSAlphabet::getMasses(size_type index) const
    {
      masses_type masses;
      
      for (const IMSElement& cit : elements_)
      {
        masses.push_back(cit.getMass(index));
      }
      return masses;
    }

    IMSAlphabet::masses_type IMSAlphabet::getAverageMasses() const
    {
      masses_type masses;
      for (const IMSElement& cit : elements_)
      {
        masses.push_back(cit.getAverageMass());
      }
      return masses;
    }

    void IMSAlphabet::sortByNames()
    {
      std::sort(elements_.begin(), elements_.end(),
                [&](const element_type& a, const element_type& b)
                { return a.getName() < b.getName(); });
    }

    void IMSAlphabet::sortByValues()
    {
      std::sort(elements_.begin(), elements_.end(), MassSortingCriteria_());
    }

    void IMSAlphabet::load(const std::string & fname)
    {
      IMSAlphabetTextParser parser;
      this->load(fname, parser);
    }

    void IMSAlphabet::load(const std::string & fname, IMSAlphabetParser<> & parser)
    {
      parser.load(fname);
      this->clear();
      for (const auto & pos : parser.getElements())
      {
        this->push_back(pos.first, pos.second);
      }
      this->sortByValues();
    }

    std::ostream & operator<<(std::ostream & os, const IMSAlphabet & alphabet)
    {
      for (IMSAlphabet::size_type i = 0; i < alphabet.size(); ++i)
      {
        os << alphabet.getElement(i) << '\n';
      }
      return os;
    }

} // namespace OpenMS // namespace ims

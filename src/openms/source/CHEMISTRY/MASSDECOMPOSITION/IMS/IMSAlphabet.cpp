// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#include <functional>
#include <algorithm>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetTextParser.h>

namespace OpenMS
{

  namespace ims
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
      const_iterator cit = elements_.begin();
      for (; cit != elements_.end(); ++cit)
      {
        if (cit->getName() == name)
        {
          return *cit;
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
      iterator it = elements_.begin();
      for (; it != elements_.end(); ++it)
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
      const_iterator cit = elements_.begin();
      for (; cit != elements_.end(); ++cit)
      {
        masses.push_back(cit->getMass(index));
      }
      return masses;
    }

    IMSAlphabet::masses_type IMSAlphabet::getAverageMasses() const
    {
      masses_type masses;
      const_iterator cit = elements_.begin();
      for (; cit != elements_.end(); ++cit)
      {
        masses.push_back(cit->getAverageMass());
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
      this->load(fname, new IMSAlphabetTextParser);
    }

    void IMSAlphabet::load(const std::string & fname, IMSAlphabetParser<> * parser)
    {
      parser->load(fname);
      this->clear();
      for (IMSAlphabetParser<>::ContainerType::const_iterator pos =
             parser->getElements().begin(),
           end = parser->getElements().end(); pos != end; ++pos)
      {
        this->push_back(pos->first, pos->second);
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

  } // namespace ims
} // namespace OpenMS

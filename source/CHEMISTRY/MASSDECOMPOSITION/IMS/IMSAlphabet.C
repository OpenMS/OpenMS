// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#include <functional>
#include <algorithm>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/compose_f_gx_t.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/compose_f_gx_hy_t.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetTextParser.h>


namespace OpenMS {

namespace ims {


const IMSAlphabet::name_type& IMSAlphabet::getName(size_type index) const
{
  return getElement(index).getName();
}


IMSAlphabet::mass_type IMSAlphabet::getMass(size_type index) const
{
  return getElement(index).getMass();
}


IMSAlphabet::mass_type IMSAlphabet::getMass(const name_type& name) const
{
  return getElement(name).getMass();
}


bool IMSAlphabet::hasName(const name_type& name) const
{
  return std::find_if(elements_.begin(), elements_.end(),
                      compose_f_gx(std::bind2nd(std::equal_to<name_type>(), name),
                                   std::mem_fun_ref(&element_type::getName))) < elements_.end();
}


const IMSAlphabet::element_type& IMSAlphabet::getElement(const name_type& name) const
{
  const_iterator cit = elements_.begin();
  for (; cit != elements_.end(); ++cit)
  {
    if (cit->getName() == name)
    {
      return *cit;
    }
  }
  throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, name + " was not found in IMSAlphabet!",String(name));
}

void IMSAlphabet::setElement(const name_type& name, mass_type mass, bool forced)
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

bool IMSAlphabet::erase(const name_type& name)
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
            compose_f_gx_hy(
              std::less<name_type>(),
              std::mem_fun_ref(&element_type::getName),
              std::mem_fun_ref(&element_type::getName)));
}


void IMSAlphabet::sortByValues()
{
  std::sort(elements_.begin(), elements_.end(), MassSortingCriteria_());
}


void IMSAlphabet::load(const std::string& fname)
{
  this->load(fname, new IMSAlphabetTextParser);
}


void IMSAlphabet::load(const std::string& fname, IMSAlphabetParser<>* parser)
{
  parser->load(fname);
  this->clear();
  for (IMSAlphabetParser<>::ContainerType::const_iterator pos =
       parser->getElements().begin(),
       end = parser->getElements().end();	pos != end; ++pos)
  {
    this->push_back(pos->first, pos->second);
  }
  this->sortByValues();
}


std::ostream& operator <<(std::ostream& os, const IMSAlphabet& alphabet)
{
  for (IMSAlphabet::size_type i = 0; i < alphabet.size(); ++i )
  {
    os << alphabet.getElement(i) << '\n';
  }
  return os;
}

} // namespace ims
} // namespace OpenMS

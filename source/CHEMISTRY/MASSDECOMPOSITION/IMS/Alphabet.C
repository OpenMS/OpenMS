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

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Alphabet.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/compose_f_gx_t.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/compose_f_gx_hy_t.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/AlphabetTextParser.h>


namespace OpenMS {

namespace ims {


const Alphabet::name_type& Alphabet::getName(size_type index) const
{
  return getElement(index).getName();
}


Alphabet::mass_type Alphabet::getMass(size_type index) const
{
  return getElement(index).getMass();
}


Alphabet::mass_type Alphabet::getMass(const name_type& name) const
{
  return getElement(name).getMass();
}


bool Alphabet::hasName(const name_type& name) const
{
  return std::find_if(elements.begin(), elements.end(),
                      compose_f_gx(std::bind2nd(std::equal_to<name_type>(), name),
                                   std::mem_fun_ref(&element_type::getName))) < elements.end();
}


const Alphabet::element_type& Alphabet::getElement(const name_type& name) const
{
  const_iterator cit = elements.begin();
  for (; cit != elements.end(); ++cit)
  {
    if (cit->getName() == name)
    {
      return *cit;
    }
  }
  throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, name + " was not found in alphabet!",String(name));
}

void Alphabet::setElement(const name_type& name, mass_type mass, bool forced)
{
  bool found = false;
  for (size_type i = 0; i < elements.size(); ++i)
  {
    if (name == elements[i].getName())
    {
      element_type element(name, mass);
      elements[i] = element;
      found = true;
      break;
    }
  }
  if (!found && forced)
  {
    this->push_back(name, mass);
  }
}

bool Alphabet::erase(const name_type& name)
{
  bool found = false;
  iterator it = elements.begin();
  for (; it != elements.end(); ++it)
  {
    if (it->getName() == name)
    {
      elements.erase(it);
      found = true;
      break;
    }
  }
  return found;
}

Alphabet::masses_type Alphabet::getMasses(size_type index) const
{
  masses_type masses;
  const_iterator cit = elements.begin();
  for (; cit != elements.end(); ++cit)
  {
    masses.push_back(cit->getMass(index));
  }
  return masses;
}


Alphabet::masses_type Alphabet::getAverageMasses() const
{
  masses_type masses;
  const_iterator cit = elements.begin();
  for (; cit != elements.end(); ++cit)
  {
    masses.push_back(cit->getAverageMass());
  }
  return masses;
}


void Alphabet::sortByNames()
{
  std::sort(elements.begin(), elements.end(),
            compose_f_gx_hy(
              std::less<name_type>(),
              std::mem_fun_ref(&element_type::getName),
              std::mem_fun_ref(&element_type::getName)));
}


void Alphabet::sortByValues()
{
  std::sort(elements.begin(), elements.end(), MassSortingCriteria());
}


void Alphabet::load(const std::string& fname)
{
  this->load(fname, new AlphabetTextParser);
}


void Alphabet::load(const std::string& fname, AlphabetParser<>* parser)
{
  parser->load(fname);
  this->clear();
  for (AlphabetParser<>::ContainerType::const_iterator pos =
       parser->getElements().begin(),
       end = parser->getElements().end();	pos != end; ++pos)
  {
    this->push_back(pos->first, pos->second);
  }
  this->sortByValues();
}


std::ostream& operator <<(std::ostream& os, const Alphabet& alphabet)
{
  for (Alphabet::size_type i = 0; i < alphabet.size(); ++i )
  {
    os << alphabet.getElement(i) << '\n';
  }
  return os;
}

} // namespace ims
} // namespace OpenMS

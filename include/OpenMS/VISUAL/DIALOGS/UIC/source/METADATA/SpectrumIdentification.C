// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SpectrumIdentification.h>

using namespace std;

namespace OpenMS {

  SpectrumIdentification::SpectrumIdentification()
    : MetaInfoInterface(),
    	id_()
  {
  }

  SpectrumIdentification::SpectrumIdentification(const SpectrumIdentification& rhs)
  	: MetaInfoInterface(rhs),
  		id_(rhs.id_),
			hits_(rhs.hits_)
  {
  }

  SpectrumIdentification::~SpectrumIdentification()
  {
  }

  SpectrumIdentification& SpectrumIdentification::operator=(const SpectrumIdentification& rhs)
  {
  	if (this == &rhs)
  	{
  		return *this;
  	}

    MetaInfoInterface::operator=(rhs);
    id_ = rhs.id_;
		hits_ = rhs.hits_;

    return *this;
  }

	// Equality operator
	bool SpectrumIdentification::operator == (const SpectrumIdentification& rhs) const
	{
		return MetaInfoInterface::operator==(rhs)
				&& id_ == rhs.id_
				&& hits_ == rhs.hits_
				;
	}

	// Inequality operator
	bool SpectrumIdentification::operator != (const SpectrumIdentification& rhs) const
	{
		return !(*this == rhs);
	}

	void SpectrumIdentification::setHits(const vector<IdentificationHit>& hits)
	{
		hits_ = hits;
	}

	void SpectrumIdentification::addHit(const IdentificationHit& hit)
	{
		hits_.push_back(hit);
	}

	const vector<IdentificationHit>& SpectrumIdentification::getHits() const
	{
		return hits_;
	}

}// namespace OpenMS

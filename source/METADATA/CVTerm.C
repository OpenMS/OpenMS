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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/CVTerm.h>

using namespace std;

namespace OpenMS
{
	// CV term implementation
	CVTerm::CVTerm()
	{
	}
			
	CVTerm::CVTerm(const String& accession, const String& name, const String& cv_identifier_ref, const String& value, const Unit& unit)
		:	accession_(accession),
			name_(name),
			cv_identifier_ref_(cv_identifier_ref),
			unit_(unit),
			value_(value)
	{
	}

	CVTerm::CVTerm(const CVTerm& rhs)
  	:	accession_(rhs.accession_),
			name_(rhs.name_),
			cv_identifier_ref_(rhs.cv_identifier_ref_),
			unit_(rhs.unit_),
			value_(rhs.value_)
	{
	}

	CVTerm::~CVTerm()
	{
	}
	
	CVTerm& CVTerm::operator = (const CVTerm& rhs)
	{
		if (this != &rhs)
		{
			accession_ = rhs.accession_;
			name_ = rhs.name_;
			cv_identifier_ref_ = rhs.cv_identifier_ref_;
			unit_ = rhs.unit_;
			value_ = rhs.value_;
		}
		return *this;
	}

	bool CVTerm::operator == (const CVTerm& rhs) const
	{
		return 	accession_ == rhs.accession_ &&
			      name_ == rhs.name_ &&
			      cv_identifier_ref_ == rhs.cv_identifier_ref_ &&
						unit_ == rhs.unit_ &&
						value_ == rhs.value_;
	}

	bool CVTerm::operator != (const CVTerm& rhs) const
	{
		return !(*this == rhs);
	}

  void CVTerm::setAccession(const String& accession)
	{
		accession_ = accession;
	}

	const String& CVTerm::getAccession() const
	{
		return accession_;
	}

  void CVTerm::setName(const String& name)
	{
		name_ = name;
	}

	const String& CVTerm::getName() const
	{
		return name_;
	}

	void CVTerm::setCVIdentifierRef(const String& cv_identifier_ref)
	{
		cv_identifier_ref_ = cv_identifier_ref;
	}

  const String& CVTerm::getCVIdentifierRef() const
	{
		return cv_identifier_ref_;
	}

	void CVTerm::setUnit(const Unit& unit)
	{
		unit_ = unit;
	}

	const CVTerm::Unit& CVTerm::getUnit() const
	{
		return unit_;
	}

	void CVTerm::setValue(const DataValue& value)
	{
		value_ = value;
	}

	const DataValue& CVTerm::getValue() const
	{
		return value_;
	}

	bool CVTerm::hasUnit() const
	{
		return unit_.accession != "";
	}

	bool CVTerm::hasValue() const
	{
		return value_ != DataValue::EMPTY;
	}

} // namespace OpenMS




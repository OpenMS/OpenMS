// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>

using namespace std;

namespace OpenMS
{
	// CV term implementation
	CVMappingTerm::CVMappingTerm()
		: use_term_name_(false),
			use_term_(false),
			is_repeatable_(false),
			allow_children_(false)
	{
	}

	CVMappingTerm::CVMappingTerm(const CVMappingTerm& rhs)
  	:	accession_(rhs.accession_),
			use_term_name_(rhs.use_term_name_),
			use_term_(rhs.use_term_),
			term_name_(rhs.term_name_),
			is_repeatable_(rhs.is_repeatable_),
			allow_children_(rhs.allow_children_),
			cv_identifier_ref_(rhs.cv_identifier_ref_)
	{
	}

	CVMappingTerm::~CVMappingTerm()
	{
	}
	
	CVMappingTerm& CVMappingTerm::operator = (const CVMappingTerm& rhs)
	{
		if (this != &rhs)
		{
			accession_ = rhs.accession_;
			use_term_name_ = rhs.use_term_name_;
			use_term_ = rhs.use_term_;
			term_name_ = rhs.term_name_;
			is_repeatable_ = rhs.is_repeatable_;
			allow_children_ =rhs.allow_children_;
			cv_identifier_ref_ = rhs.cv_identifier_ref_;
		}
		return *this;
	}

	bool CVMappingTerm::operator == (const CVMappingTerm& rhs) const
	{
		return 	accession_ == rhs.accession_ &&
			      use_term_name_ == rhs.use_term_name_ &&
			      use_term_ == rhs.use_term_ &&
			      term_name_ == rhs.term_name_ &&
						is_repeatable_ == rhs.is_repeatable_ &&
			      allow_children_ == rhs.allow_children_ &&
			      cv_identifier_ref_ == rhs.cv_identifier_ref_;
	}

	bool CVMappingTerm::operator != (const CVMappingTerm& rhs) const
	{
		return !(*this == rhs);
	}

  void CVMappingTerm::setAccession(const String& accession)
	{
		accession_ = accession;
	}

	const String& CVMappingTerm::getAccession() const
	{
		return accession_;
	}

  void CVMappingTerm::setUseTermName(bool use_term_name)
	{
		use_term_name_ = use_term_name;
	}

  bool CVMappingTerm::getUseTermName() const
	{
		return use_term_name_;
	}

  void CVMappingTerm::setUseTerm(bool use_term)
	{
		use_term_ = use_term;
	}

	bool CVMappingTerm::getUseTerm() const
	{
		return use_term_;
	}

  void CVMappingTerm::setTermName(const String& term_name)
	{
		term_name_ = term_name;
	}

	const String& CVMappingTerm::getTermName() const
	{
		return term_name_;
	}

  void CVMappingTerm::setIsRepeatable(bool is_repeatable)
	{
		is_repeatable_ = is_repeatable;
	}

  bool CVMappingTerm::getIsRepeatable() const
	{
		return is_repeatable_;
	}

  void CVMappingTerm::setAllowChildren(bool allow_children)
	{
		allow_children_ = allow_children;
	}

  bool CVMappingTerm::getAllowChildren() const
	{
		return allow_children_;
	}

	void CVMappingTerm::setCVIdentifierRef(const String& cv_identifier_ref)
	{
		cv_identifier_ref_ = cv_identifier_ref;
	}

  const String& CVMappingTerm::getCVIdentifierRef() const
	{
		return cv_identifier_ref_;
	}

} // namespace OpenMS




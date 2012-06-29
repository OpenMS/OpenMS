// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/CVTermList.h>

using namespace std;

namespace OpenMS
{
	// CV term implementation
	CVTermList::CVTermList()
		:	MetaInfoInterface()
	{
	}

	CVTermList::CVTermList(const CVTermList& rhs)
  	:	MetaInfoInterface(rhs),
			cv_terms_(rhs.cv_terms_)
	{
	}

	CVTermList::~CVTermList()
	{
	}
	
	CVTermList& CVTermList::operator = (const CVTermList& rhs)
	{
		if (this != &rhs)
		{
			MetaInfoInterface::operator = (rhs);
			cv_terms_ = rhs.cv_terms_;
		}
		return *this;
	}

  void CVTermList::addCVTerm(const CVTerm& cv_term)
	{
		// TODO exception if empty
		cv_terms_[cv_term.getAccession()].push_back(cv_term);
	}

	void CVTermList::setCVTerms(const vector<CVTerm>& cv_terms)
	{
		for (vector<CVTerm>::const_iterator it = cv_terms.begin(); it != cv_terms.end(); ++it)
		{
			addCVTerm(*it);
		}
		return;
	}

	void CVTermList::replaceCVTerm(const CVTerm& cv_term)
  {
    std::vector<CVTerm> tmp;
    tmp.push_back(cv_term);
    cv_terms_[cv_term.getAccession()] = tmp;
  }

	void CVTermList::replaceCVTerms(const vector<CVTerm>& cv_terms, const String& accession)
	{
		cv_terms_[accession] = cv_terms;
	}

	void CVTermList::replaceCVTerms(const Map<String, vector<CVTerm> >& cv_term_map)
	{
		cv_terms_ = cv_term_map;
	}

	const Map<String, vector<CVTerm> >& CVTermList::getCVTerms() const
	{
		return cv_terms_;
	}

	bool CVTermList::hasCVTerm(const String& accession) const
	{
		return cv_terms_.has(accession);
	}

	bool CVTermList::operator == (const CVTermList& cv_term_list) const
	{
		return MetaInfoInterface::operator == (cv_term_list) && cv_terms_ == cv_term_list.cv_terms_;
	}

	bool CVTermList::operator != (const CVTermList& cv_term_list) const
	{
		return !(*this == cv_term_list);
	}

	bool CVTermList::empty() const
	{
		return cv_terms_.empty();
	}

} // namespace OpenMS




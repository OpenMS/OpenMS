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

#include <OpenMS/DATASTRUCTURES/CVMappings.h>

using namespace std;

namespace OpenMS
{
	CVMappings::CVMappings()
	{
	}

	CVMappings::CVMappings(const CVMappings& rhs)
		: mapping_rules_(rhs.mapping_rules_),
			cv_references_(rhs.cv_references_),
			cv_references_vector_(rhs.cv_references_vector_)
	{
	}

	CVMappings::~CVMappings()
	{
	}
	
	CVMappings& CVMappings::operator = (const CVMappings& rhs)
	{
		if (this != &rhs)
		{
			mapping_rules_ = rhs.mapping_rules_;
			cv_references_ = rhs.cv_references_;
			cv_references_vector_ = rhs.cv_references_vector_;
		}
		return *this;
	}

	bool CVMappings::operator == (const CVMappings& rhs) const
	{
		return 	mapping_rules_ == rhs.mapping_rules_ && 
						cv_references_ == rhs.cv_references_ &&
						cv_references_vector_ == rhs.cv_references_vector_;
	}

	bool CVMappings::operator != (const CVMappings& rhs) const
	{
		return !(*this == rhs);
	}

  void CVMappings::setMappingRules(const vector<CVMappingRule>& cv_mapping_rules)
	{
		mapping_rules_ = cv_mapping_rules;
	}

  const vector<CVMappingRule>& CVMappings::getMappingRules() const
	{
		return mapping_rules_;
	}

  void CVMappings::addMappingRule(const CVMappingRule& cv_mapping_rule)
	{
		mapping_rules_.push_back(cv_mapping_rule);
	}

  void CVMappings::setCVReferences(const vector<CVReference>& cv_references)
	{
		for (vector<CVReference>::const_iterator it = cv_references.begin(); it != cv_references.end(); ++it)
		{
			cv_references_[it->getIdentifier()] = *it;
			cv_references_vector_.push_back(*it);
		}
	}

	const vector<CVReference>& CVMappings::getCVReferences() const
	{
		return cv_references_vector_;
	}

  void CVMappings::addCVReference(const CVReference& cv_reference)
	{
		if (hasCVReference(cv_reference.getIdentifier()))
		{
			cerr << "CVMappings: Warning: CV reference with identifier '" << cv_reference.getIdentifier() << "' already existing, ignoring it!" << endl;
			return;
		}
		cv_references_[cv_reference.getIdentifier()] = cv_reference;
		cv_references_vector_.push_back(cv_reference);
	}

  bool CVMappings::hasCVReference(const String& identifier)
	{
		return cv_references_.has(identifier);
	}
	
} // namespace OpenMS




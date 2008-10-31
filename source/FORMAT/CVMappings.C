// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CVMappings.h>

using namespace std;

namespace OpenMS
{
	// CV term implementation
	CVMappings::CVTerm::CVTerm()
		: use_term_name_(false),
			use_term_(false),
			is_repeatable_(false),
			allow_children_(false)
	{
	}

	CVMappings::CVTerm::CVTerm(const CVTerm& rhs)
  	:	accession_(rhs.accession_),
			use_term_name_(rhs.use_term_name_),
			use_term_(rhs.use_term_),
			term_name_(rhs.term_name_),
			is_repeatable_(rhs.is_repeatable_),
			allow_children_(rhs.allow_children_),
			cv_identifier_ref_(rhs.cv_identifier_ref_)
	{
	}

	CVMappings::CVTerm::~CVTerm()
	{
	}
	
	CVMappings::CVTerm& CVMappings::CVTerm::operator = (const CVTerm& rhs)
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

  void CVMappings::CVTerm::setAccession(const String& accession)
	{
		accession_ = accession;
	}

	const String& CVMappings::CVTerm::getAccession() const
	{
		return accession_;
	}

  void CVMappings::CVTerm::setUseTermName(bool use_term_name)
	{
		use_term_name_ = use_term_name;
	}

  bool CVMappings::CVTerm::getUseTermName() const
	{
		return use_term_name_;
	}

  void CVMappings::CVTerm::setUseTerm(bool use_term)
	{
		use_term_ = use_term;
	}

	bool CVMappings::CVTerm::getUseTerm() const
	{
		return use_term_;
	}

  void CVMappings::CVTerm::setTermName(const String& term_name)
	{
		term_name_ = term_name;
	}

	const String& CVMappings::CVTerm::getTermName() const
	{
		return term_name_;
	}

  void CVMappings::CVTerm::setIsRepeatable(bool is_repeatable)
	{
		is_repeatable_ = is_repeatable;
	}

  bool CVMappings::CVTerm::getIsRepeatable() const
	{
		return is_repeatable_;
	}

  void CVMappings::CVTerm::setAllowChildren(bool allow_children)
	{
		allow_children_ = allow_children;
	}

  bool CVMappings::CVTerm::getAllowChildren() const
	{
		return allow_children_;
	}

	void CVMappings::CVTerm::setCVIdentifierRef(const String& cv_identifier_ref)
	{
		cv_identifier_ref_ = cv_identifier_ref;
	}

  const String& CVMappings::CVTerm::getCVIdentifierRef() const
	{
		return cv_identifier_ref_;
	}


	// CV mapping rule implementation
	CVMappings::CVMappingRule::CVMappingRule()
		: requirement_level_(CVMappingRule::MUST),
			combinations_logic_(CVMappingRule::OR)
	{
	}

	CVMappings::CVMappingRule::CVMappingRule(const CVMappingRule& rhs)
		:	identifier_(rhs.identifier_),
			element_path_(rhs.element_path_),
      requirement_level_(rhs.requirement_level_),
      scope_path_(rhs.scope_path_),
      combinations_logic_(rhs.combinations_logic_),
      cv_terms_(rhs.cv_terms_)
	{
	}

	CVMappings::CVMappingRule::~CVMappingRule()
	{
	}
	
	CVMappings::CVMappingRule& CVMappings::CVMappingRule::operator = (const CVMappingRule& rhs)
	{
		if (this != &rhs)
		{
			identifier_ = rhs.identifier_;
      element_path_ = rhs.element_path_;
      requirement_level_ = rhs.requirement_level_;
      scope_path_ = rhs.scope_path_;
      combinations_logic_ = rhs.combinations_logic_;
      cv_terms_ = rhs.cv_terms_;
		}
		return *this;
	}

  void CVMappings::CVMappingRule::setIdentifier(const String& identifier)
	{
		identifier_ = identifier;
	}

	const String& CVMappings::CVMappingRule::getIdentifier() const
	{
		return identifier_;
	}

	void CVMappings::CVMappingRule::setElementPath(const String& element_path)
	{
		element_path_ = element_path;
	}

	const String& CVMappings::CVMappingRule::getElementPath() const
	{
		return element_path_;
	}

	void CVMappings::CVMappingRule::setRequirementLevel(RequirementLevel level)
	{
		requirement_level_ = level;
	}

	CVMappings::CVMappingRule::RequirementLevel CVMappings::CVMappingRule::getRequirementLevel() const
	{
		return requirement_level_;
	}

	void CVMappings::CVMappingRule::setCombinationsLogic(CombinationsLogic combinations_logic)
	{
		combinations_logic_ = combinations_logic;
	}
	
	CVMappings::CVMappingRule::CombinationsLogic CVMappings::CVMappingRule::getCombinationsLogic() const
	{
		return combinations_logic_;
	}

	void CVMappings::CVMappingRule::setScopePath(const String& path)
	{
		scope_path_ = path;
	}

	const String& CVMappings::CVMappingRule::getScopePath() const
	{
		return scope_path_;
	}

	void CVMappings::CVMappingRule::setCVTerms(const vector<CVTerm>& cv_terms)
	{
		cv_terms_ = cv_terms;
	}
	 
	const vector<CVMappings::CVTerm>& CVMappings::CVMappingRule::getCVTerms() const
	{
		return cv_terms_;
	}
	
  void CVMappings::CVMappingRule::addCVTerm(const CVMappings::CVTerm& cv_term)
	{
		cv_terms_.push_back(cv_term);
	}

	// CV reference implementation
	CVMappings::CVReference::CVReference()
	{
	}

	CVMappings::CVReference::~CVReference()
	{
	}
	
	CVMappings::CVReference::CVReference(const CVReference& rhs)
		: name_(rhs.name_),
			identifier_(rhs.identifier_)
	{
	}

	CVMappings::CVReference& CVMappings::CVReference::operator = (const CVReference& rhs)
	{
		if (this != &rhs)
		{
			name_ = rhs.name_;
			identifier_ = rhs.identifier_;
		}
		return *this;
	}

	void CVMappings::CVReference::setName(const String& name)
	{
		name_ = name;
	}

	const String& CVMappings::CVReference::getName() const
	{
		return name_;
	}

	void CVMappings::CVReference::setIdentifier(const String& identifier)
	{
		identifier_ = identifier;
	}

	const String& CVMappings::CVReference::getIdentifier() const
	{
		return identifier_;
	}

	// CV mappings implementation
	CVMappings::CVMappings()
	{
	}

	CVMappings::CVMappings(const CVMappings& rhs)
		: mapping_rules_(rhs.mapping_rules_),
			cv_references_(rhs.cv_references_)
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
		}
		return *this;
	}

  void CVMappings::setMappingRules(const vector<CVMappingRule>& cv_mapping_rules)
	{
		mapping_rules_ = cv_mapping_rules;
	}

  const vector<CVMappings::CVMappingRule>& CVMappings::getMappingRules() const
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
		}
	}

  void CVMappings::addCVReference(const CVReference& cv_reference)
	{
		if (hasCVReference(cv_reference.getIdentifier()))
		{
			cerr << "CVMappings: Warning: CV reference with identifier '" << cv_reference.getIdentifier() << "' already existing, ignoreing it!" << endl;
			return;
		}
		cv_references_[cv_reference.getIdentifier()] = cv_reference;
	}

  bool CVMappings::hasCVReference(const String& identifier)
	{
		return cv_references_.has(identifier);
	}
	
} // namespace OpenMS




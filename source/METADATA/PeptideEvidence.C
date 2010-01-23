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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideEvidence.h>

#include <algorithm>

using namespace std;

namespace OpenMS 
{
	// default constructor
  PeptideEvidence::PeptideEvidence()
		:	MetaInfoInterface(),
      db_sequence_ref_(""),
      translation_table_ref_(""),
      start_(-1),
      end_(-1),
      pre_('?'),
      post_('?'),
      id_(""),
      name_(""),
      missed_cleavages_(-1),
      is_decoy_(false),
			frame_(-1)
  {
  }
  
	// copy constructor
  PeptideEvidence::PeptideEvidence(const PeptideEvidence& rhs)
		:	MetaInfoInterface(rhs),
      db_sequence_ref_(rhs.db_sequence_ref_),
      translation_table_ref_(rhs.translation_table_ref_),
      start_(rhs.start_),
      end_(rhs.end_),
      pre_(rhs.pre_),
      post_(rhs.post_),
      id_(rhs.id_),
      name_(rhs.name_),
      missed_cleavages_(rhs.missed_cleavages_),
      is_decoy_(rhs.is_decoy_),
      frame_(rhs.frame_)
  {
  }
  
	// destructor
  PeptideEvidence::~PeptideEvidence()
  {  	
  }
   
  PeptideEvidence& PeptideEvidence::operator= (const PeptideEvidence& rhs)
  {
	 	if (this == &rhs)
  	{
  		return *this;
  	}
  	
  	MetaInfoInterface::operator=(rhs);
    db_sequence_ref_ = rhs.db_sequence_ref_;
    translation_table_ref_ = rhs.translation_table_ref_;
    start_ = rhs.start_;
    end_ = rhs.end_;
    pre_ = rhs.pre_;
    post_ = rhs.post_;
    id_ = rhs.id_;
    name_ = name_;
    missed_cleavages_ = rhs.missed_cleavages_;
    is_decoy_ = rhs.is_decoy_;
    frame_ = rhs.frame_;

    return *this;
  }
	
	bool PeptideEvidence::operator == (const PeptideEvidence& rhs) const	
	{
		return MetaInfoInterface::operator==(rhs) &&
	    db_sequence_ref_ == rhs.db_sequence_ref_ &&
  	  translation_table_ref_ == rhs.translation_table_ref_ &&
  	  start_ == rhs.start_ &&
  	  end_ == rhs.end_ &&
  	  pre_ == rhs.pre_ &&
    	post_ == rhs.post_ &&
   	 	id_ == rhs.id_ &&
    	name_ == name_ &&
    	missed_cleavages_ == rhs.missed_cleavages_ &&
    	is_decoy_ == rhs.is_decoy_ &&
    	frame_ == rhs.frame_;
	}

	bool PeptideEvidence::operator != (const PeptideEvidence& rhs) const	
	{
		return !operator==(rhs);
	}

      
	const String& PeptideEvidence::getDBSequenceRef() const
	{
		return db_sequence_ref_;
	}
	
	void PeptideEvidence::setDBSequenceRef(const String& rhs)
	{
		db_sequence_ref_ = rhs;
	}

	const String& PeptideEvidence::getTranslationTableRef() const
	{
		return translation_table_ref_;
	}

	void PeptideEvidence::setTranslationTableRef(const String&  rhs)
	{
		translation_table_ref_ = rhs;
	}


	void PeptideEvidence::setStart(Int start)
	{
		start_ = start;
	}

	Int PeptideEvidence::getStart() const
	{
		return start_;
	}

 	void PeptideEvidence::setEnd(Int end)
	{
		end_ = end;
	}

	Int PeptideEvidence::getEnd() const
	{
		return end_;
	}

	void PeptideEvidence::setPre(char rhs)
	{
		pre_ = rhs;
	}

  char PeptideEvidence::getPre() const
	{
		return pre_;
	}

	void PeptideEvidence::setPost(char rhs)
	{
		post_ = rhs;
	}

	char PeptideEvidence::getPost() const
	{
		return post_;
	}

	void PeptideEvidence::setId(const String& id)
	{
		id_ = id;
	}
	
	const String& PeptideEvidence::getId() const
	{
		return id_;
	}

	void PeptideEvidence::setName(const String& name)
	{
		name_ = name;
	}
      
	const String& PeptideEvidence::getName() const
	{
		return name_;
	}

  void PeptideEvidence::setMissedCleavages(Int rhs)
	{
		missed_cleavages_ = rhs;
	}

  Int PeptideEvidence::getMissedCleavages() const
	{
		return missed_cleavages_;
	}

	void PeptideEvidence::setIsDecoy(bool is_decoy)
	{
		is_decoy_ = is_decoy;
	}

  bool PeptideEvidence::getIsDecoy() const
	{
		return is_decoy_;
	}

	void PeptideEvidence::setFrame(Int frame)
	{
		frame_ = frame;
	}

	Int PeptideEvidence::getFrame() const
	{
		return frame_;
	}


} // namespace OpenMS

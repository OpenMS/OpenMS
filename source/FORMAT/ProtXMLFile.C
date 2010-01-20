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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ProtXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <fstream>
#include <limits>

using namespace std;

namespace OpenMS 
{

	ProtXMLFile::ProtXMLFile()
		: XMLHandler("","1.2"),
			XMLFile("/SCHEMAS/protXML_v6.xsd","6.0")
	{
	}

  void ProtXMLFile::load(const String& filename,  ProteinIdentification& protein_ids, PeptideIdentification& peptide_ids)
  {
  	//Filename for error messages in XMLHandler
  	file_ = filename;

		resetMembers_();
  	
  	// reset incoming data
  	protein_ids = ProteinIdentification();
  	peptide_ids = PeptideIdentification();
  	
  	// remember data link while parsing
  	prot_id_ = &protein_ids;
  	pep_id_ = &peptide_ids;
		
		parse_(filename,this);
  }
  					 
  void ProtXMLFile::store(const String& /*filename*/, const ProteinIdentification& /*protein_ids*/, const PeptideIdentification& /*peptide_ids*/, const String& /*document_id*/)
  {
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		// resetMembers_();
  }


	/// reset members
	void ProtXMLFile::resetMembers_()
	{
		prot_id_ = 0;
		pep_id_ = 0;
		pep_hit_ = 0;
		protein_name_to_index_.clear();
		protein_group_ = ProteinGroup();
		protein_tag_count_ = 0;
		master_protein_index_=0;
	}
  
	void ProtXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{		
		String tag = sm_.convert(qname);
		
		// identifier for Protein & PeptideIdentification
		//<program_details analysis="proteinprophet" time="2009-11-29T18:30:03"
		
		if (tag =="program_details")
		{
			String id = attributeAsString_(attributes,"analysis");
			id += "_" + String(attributeAsString_(attributes,"time"));
			
			prot_id_->setIdentifier(id);
			pep_id_->setIdentifier(id);
		}

		if (tag =="protein_group")
		{
			protein_group_ = ProteinGroup();
			protein_group_.id = attributeAsString_(attributes,"group_number");
			protein_group_.probability = attributeAsDouble_(attributes,"probability");
			protein_tag_count_ = 0;
		}
		else if (tag =="protein")
		{
			++protein_tag_count_;
			if (protein_tag_count_>1)
			{ // We currently assume that all proteins in a protein_group are indistinguisable.
			  // A second protein-tag indicates that this is not the case. We thus might need to introduce subgroups of indistinguishable proteins within a group
				LOG_ERROR << "ProtXMLFile parsing: unexpected second <protein> tag. Internal protein group structure might need addtional information to represent this!\n";
			}
			
			// see if the Protein is already known (just a precaution)
			String protein_name = attributeAsString_(attributes,"protein_name");
			registerProtein_(protein_name); // create new protein (if required)
			
			// remember last real protein
			master_protein_index_ = protein_name_to_index_[protein_name];

			// fill protein with life
			prot_id_->getHits()[master_protein_index_].setCoverage(attributeAsDouble_(attributes,"percent_coverage"));
			prot_id_->getHits()[master_protein_index_].setScore(attributeAsDouble_(attributes,"probability"));
			
		}
		else if (tag =="indistinguishable_protein")
		{
			// see if the Protein is already known (just a precaution)
			String protein_name = attributeAsString_(attributes,"protein_name");
			registerProtein_(protein_name);
			// score of group-leader might not be transferrable (due to protein-length etc),
			// so we set it to -1
			prot_id_->getHits()[protein_name_to_index_[protein_name]].setScore(-1);
			
		}
		else if (tag =="peptide")
		{
			// If a peptide is degenerate it will show in multiple groups, but have different statistics (e.g. 'nsp_adjusted_probability')
			// We thus treat each instance as a separate peptide
			// todo/improvement: link them by a group in PeptideIdentification?!
			
			pep_hit_ = new PeptideHit;
			
			pep_hit_->setSequence(attributeAsString_(attributes,"peptide_sequence"));
			pep_hit_->setCharge(attributeAsInt_(attributes,"charge"));
			pep_hit_->setScore(attributeAsDouble_(attributes,"nsp_adjusted_probability"));
			// even if the peptide is degenerate, we only store the protein it is found at 
			// (for the other proteins, the peptides' attributes will be different)
			pep_hit_->addProteinAccession(prot_id_->getHits()[master_protein_index_].getAccession());
			pep_hit_->setMetaValue("is_unique", String(attributeAsString_(attributes,"is_nondegenerate_evidence"))=="Y" ? 1 : 0);
			pep_hit_->setMetaValue("is_contributing", String(attributeAsString_(attributes,"is_contributing_evidence"))=="Y" ? 1 : 0);
		}
		else if (tag =="mod_aminoacid_mass")
		{ 
			// relates to the last seen peptide (we hope)
			Size position = attributeAsInt_(attributes,"position");
			DoubleReal mass = attributeAsDouble_(attributes,"mass");
			AASequence temp_aa_sequence(pep_hit_->getSequence());

			String temp_description = "";
			matchModification_(mass, temp_aa_sequence[position - 1].getOneLetterCode() , temp_description);

			// e.g. Carboxymethyl (C)
			vector<String> mod_split;
			temp_description.split(' ', mod_split);
			if (mod_split.size() == 2)
			{
				if (mod_split[1] == "(C-term)" || ModificationsDB::getInstance()->getModification(temp_description).getTermSpecificity() == ResidueModification::C_TERM)
				{
					temp_aa_sequence.setCTerminalModification(mod_split[0]);
				}
				else
				{
					if (mod_split[1] == "(N-term)" || ModificationsDB::getInstance()->getModification(temp_description).getTermSpecificity() == ResidueModification::N_TERM)
					{
						temp_aa_sequence.setNTerminalModification(mod_split[0]);
					}
					else
					{
						// search this mod, if not directly use a general one
						temp_aa_sequence.setModification(position - 1, mod_split[0]);
					}
				}
			}
			else
			{
				error(LOAD, String("Cannot parse modification '") + temp_description + "@" + position + "'");
			}

			pep_hit_->setSequence(temp_aa_sequence);
		}
	}
	
	void ProtXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		String tag = sm_.convert(qname);
		
		
		if (tag == "protein_group")
		{
		  prot_id_->insertGroup(protein_group_);
		}
		else if (tag =="peptide")
		{
			pep_id_->insertHit(*pep_hit_);
			delete pep_hit_;
		}
	}
	
	void ProtXMLFile::registerProtein_(const String& protein_name)
	{
			if (protein_name_to_index_.has(protein_name))
			{// existing protein (dangerous)
				std::cout << "ProtXMLFile parsing: protein '" << protein_name << "' already seen elsewhere in this document. Skipping creation, just referencing it. This might result in overwritten attributes!\n";
			}
			else
			{// new protein
				prot_id_->insertHit(ProteinHit());
				protein_name_to_index_[protein_name] = prot_id_->getHits().size()-1;
				prot_id_->getHits().back().setAccession(protein_name);
			}
			
			// Add protein to group
			protein_group_.indices.insert(protein_name_to_index_[protein_name]);

	}
	
	void ProtXMLFile::matchModification_(const DoubleReal mass, const String& origin, String& modification_description)
	{
		DoubleReal mod_mass = mass - ResidueDB::getInstance()->getResidue(origin)->getMonoWeight(Residue::Internal);
		vector<String> mods;
		ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, origin, mod_mass, 0.001);
			
		if (mods.size() == 1)
    {
			modification_description = mods[0];
    }
    else
    {
     	if (mods.size() == 0)
      {
      	error(LOAD, String("Cannot find modification '") + String(mass) + " " + String(origin) + "'");
      }
      else
      {
       	String mod_str = mods[0];
        for (vector<String>::const_iterator mit = ++mods.begin(); mit != mods.end(); ++mit)
        {
        	mod_str += ", " + *mit;
        }
				error(LOAD, "Modification '" + String(mass) + "' is not uniquely defined by the given data. Using '" + mods[0] +  "' to represent any of '" + mod_str + "'!");
				modification_description = mods[0];
      }
		}
	}

} // namespace OpenMS

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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

using namespace std;

namespace OpenMS 
{

	OMSSAXMLFile::OMSSAXMLFile()
		:	XMLHandler("", 1.1),
			XMLFile(),
			peptide_identifications_(0)
	{
		readMappingFile_();
	}
	
	OMSSAXMLFile::~OMSSAXMLFile()
	{
	}
	
  void OMSSAXMLFile::load(const String& filename, ProteinIdentification& protein_identification, vector<PeptideIdentification>& peptide_identifications, bool load_proteins)
  {		
		//Filename for error messages in XMLHandler
    file_ = filename;

		load_proteins_ = load_proteins;
		peptide_identifications_ = &peptide_identifications;

		parse_(filename, this);


		DateTime now = DateTime::now();
		String identifier("OMSSA_" + now.get());
	
		// post-processing
		vector<String> accessions;
		for (vector<PeptideIdentification>::iterator it = peptide_identifications.begin(); it != peptide_identifications.end(); ++it)
		{
			it->setScoreType("OMSSA");
			it->setHigherScoreBetter(false);
			it->setIdentifier(identifier);
			it->assignRanks();

			if (load_proteins)
			{
				for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
				{
					accessions.insert(accessions.end(), pit->getProteinAccessions().begin(), pit->getProteinAccessions().end());
				}
			}
		}

		if (load_proteins)
		{
			sort(accessions.begin(), accessions.end());
			vector<String>::const_iterator end_unique = unique(accessions.begin(), accessions.end());

			for (vector<String>::const_iterator it = accessions.begin(); it != end_unique; ++it)
			{
				ProteinHit hit;
				hit.setAccession(*it);
				protein_identification.insertHit(hit);
			}
		

			// E-values
			protein_identification.setHigherScoreBetter(false);
			protein_identification.setScoreType("OMSSA");
			protein_identification.setIdentifier(identifier);
		}
		
		// version of OMSSA is not available
		// Date of the search is not available -> set it to now
		protein_identification.setDateTime(now);
		protein_identification.setIdentifier(identifier);

		// search parameters are also not available
  }  					 

  void OMSSAXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& /*attributes*/)
	{
		tag_ = String(sm_.convert(qname)).trim();

		if (tag_ == "MSHitSet_number")
		{
			return;
		}

		if (tag_ == "MSPepHit")
		{
			return;
		}
		if (tag_ == "MSHits")
		{
			return;
		}

		if (tag_ == "MSHitSet")
		{
			return;
		}
	}
	  
  void OMSSAXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
 	{
 		tag_ = String(sm_.convert(qname)).trim();

		// protein hits (MSPepHits) are handled in ::characters(...)

		// end of peptide hit 
 		if (tag_ == "MSHits")
		{
			actual_peptide_id_.insertHit(actual_peptide_hit_);
			actual_peptide_hit_ = PeptideHit();
			tag_ = "";
			return;
		}

		// end of peptide id
		if (tag_ == "MSHitSet")
		{
			peptide_identifications_->push_back(actual_peptide_id_);
			actual_peptide_id_.assignRanks();
			actual_peptide_id_ = PeptideIdentification();
		}

		/*
		Modifications:
		
 		<MSModHit>
			<MSModHit_site>1</MSModHit_site>
      <MSModHit_modtype>
      	<MSMod>3</MSMod>
      </MSModHit_modtype>
    </MSModHit>	  
		*/
		if (tag_ == "MSModHit")
		{
			if (mods_map_.has(actual_mod_type_.toInt()) && mods_map_[actual_mod_type_.toInt()].size() > 0)
			{
				if (mods_map_[actual_mod_type_.toInt()].size() > 1)
				{
					warning(STORE, String("Cannot determine exact type of modification of position ") + actual_mod_site_ + " in sequence " + actual_peptide_hit_.getSequence().toString() + " using modification " + actual_mod_type_ + " - using first possibility!");
				}
				AASequence pep = actual_peptide_hit_.getSequence();
				pep.setModification(actual_mod_site_, mods_map_[actual_mod_type_.toInt()].begin()->getFullName());
				actual_peptide_hit_.setSequence(pep);
			}
			else
			{
				warning(STORE, String("Cannot find PSI-MOD mapping for mod -  ingoring '") + actual_mod_type_ + "'");
			}
		}
		
		tag_ = "";
 	} 

  void OMSSAXMLFile::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
  {
		String value = ((String)sm_.convert(chars)).trim();
		// MSPepHit section
    // <MSPepHit_start>0</MSPepHit_start>
    // <MSPepHit_stop>8</MSPepHit_stop>
    // <MSPepHit_accession>6599</MSPepHit_accession>
    // <MSPepHit_defline>CRHU2 carbonate dehydratase (EC 4.2.1.1) II [validated] - human</MSPepHit_defline>
    // <MSPepHit_protlength>260</MSPepHit_protlength>
    // <MSPepHit_oid>6599</MSPepHit_oid>
		if (tag_ == "MSPepHit_start")
		{
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_stop")
		{
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_accession")
		{
			if (load_proteins_)
			{
      	actual_peptide_hit_.addProteinAccession(value);
			}
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_defline")
		{
			// TODO add defline to ProteinHit?
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_protlength")
		{
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_oid")
		{
      tag_ = "";
      return;
		}

		// MSHits section
		// <MSHits_evalue>0.00336753988893542</MSHits_evalue>
		// <MSHits_pvalue>1.30819399070598e-08</MSHits_pvalue>
		// <MSHits_charge>1</MSHits_charge>
		// <MSHits_pepstring>MSHHWGYGK</MSHits_pepstring>
    // <MSHits_mass>1101492</MSHits_mass>
    // <MSHits_pepstart></MSHits_pepstart>
    // <MSHits_pepstop>H</MSHits_pepstop>
    // <MSHits_theomass>1101484</MSHits_theomass>
		if (tag_ == "MSHits_evalue")
		{
			actual_peptide_hit_.setScore(value.toDouble());
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_charge")
		{
			actual_peptide_hit_.setCharge(value.toInt());
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pvalue")
		{
			// TODO extra field?
			//actual_peptide_hit_.setScore(value.toDouble());
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pepstring")
		{
			AASequence seq = value.trim();
			if (mod_def_set_.getNumberOfFixedModifications() != 0 && seq.isValid())
			{
				set<String> fixed_mod_names = mod_def_set_.getFixedModificationNames();
				for (set<String>::const_iterator it = fixed_mod_names.begin(); it != fixed_mod_names.end(); ++it)
				{
					String origin = ModificationsDB::getInstance()->getModification(*it).getOrigin();
					UInt position(0);
					for (AASequence::Iterator ait = seq.begin(); ait != seq.end(); ++ait, ++position)
					{
						if (ait->getOneLetterCode() == origin)
						{
							seq.setModification(position, *it);
						}
					}
				}
			}
			actual_peptide_hit_.setSequence(seq);
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_mass")
		{
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pepstart")
		{
			if (value != "")
			{
				actual_peptide_hit_.setAABefore(value[0]);
			}
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pepstop")
		{
			if (value != "")
			{
				actual_peptide_hit_.setAAAfter(value[0]);
			}
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_theomass")
		{
			
			tag_ = "";
			return;
		}

		// modifications
		//<MSHits_mods>
    //	<MSModHit>
		//		<MSModHit_site>4</MSModHit_site>
		//		<MSModHit_modtype>
		//			<MSMod>2</MSMod>
		//		</MSModHit_modtype>
		//	</MSModHit>
		//</MSHits_mods>

		
		if (tag_ == "MSHits_mods")
		{
			actual_mod_site_ = 0;
			actual_mod_type_ = "";
		}
		if (tag_ == "MSModHit_site")
		{
			actual_mod_site_ = value.trim().toInt();
		}
		if (tag_ == "MSMod")
		{
			actual_mod_type_ = value.trim();
		}
	
		// m/z value and rt 
		if (tag_ == "MSHitSet_ids_E")
		{
			if (value.trim() != "")
			{
				if (value.has('_'))
				{
					String mz(value.prefix('_'));
					String rt(value.suffix('_'));
					try
					{
						actual_peptide_id_.setMetaValue("MZ", mz.toDouble());
						actual_peptide_id_.setMetaValue("RT", rt.toDouble());				
					}
					catch (...)
					{
						// if exception happens to occur here, s.th. went wrong, e.g. the value does not contains numbers
					}
				}
			}
		}
	}

	void OMSSAXMLFile::readMappingFile_()
	{
		String file = File::find("CHEMISTRY/OMSSA_modification_mapping");
		TextFile infile(file);
		
		for (TextFile::ConstIterator it = infile.begin(); it != infile.end(); ++it)
		{
			vector<String> split;
			it->split(',', split);

			if (it->size() > 0 && (*it)[0] != '#')
			{
				if (split.size() < 2)
				{
					fatalError(LOAD, String("Invalid mapping file line: '") + *it + "'");
				}
				vector<ResidueModification> mods;
				for (Size i = 2; i != split.size(); ++i)
				{
					String tmp(split[i].trim());
					if (tmp.size() != 0)
					{
						mods.push_back(ModificationsDB::getInstance()->getModification(tmp));
					}
				}
				mods_map_[split[0].trim().toInt()] = mods;
			}
		}
	}

	void OMSSAXMLFile::setModificationDefinitionsSet(const ModificationDefinitionsSet& mod_set)
	{
		mod_def_set_ = mod_set;
	}
	
} // namespace OpenMS

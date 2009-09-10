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

#include <OpenMS/FORMAT/HANDLERS/MzIdentMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <set>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{

  MzIdentMLHandler::MzIdentMLHandler(const Identification& id, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			id_(0),
			cid_(&id)
  {
  	cv_.loadFromOBO("PI",File::find("/CV/psi-ms.obo"));
  }

  MzIdentMLHandler::MzIdentMLHandler(Identification& id, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			id_(&id),
			cid_(0)
  {
  	cv_.loadFromOBO("MS",File::find("/CV/psi-ms.obo"));
  }	

	MzIdentMLHandler::~MzIdentMLHandler()
	{
	}

	void MzIdentMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		tag_ = sm_.convert(qname);

		open_tags_.push_back(tag_);

		//determine parent tag
    String parent_tag;
    if (open_tags_.size() > 1)
		{
			 parent_tag = *(open_tags_.end()-2);
		}
    String parent_parent_tag;
    if (open_tags_.size() > 2)
		{
			parent_parent_tag = *(open_tags_.end()-3);
		}

		static const XMLCh* s_value = xercesc::XMLString::transcode("value");
    static const XMLCh* s_unit_accession = xercesc::XMLString::transcode("unitAccession");
    static const XMLCh* s_cv_ref = xercesc::XMLString::transcode("cvRef");
    static const XMLCh* s_name = xercesc::XMLString::transcode("name");
    static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");


		if (tag_ == "cvParam")
		{
      String value, unit_accession, cv_ref;
      optionalAttributeAsString_(value, attributes, s_value);
      optionalAttributeAsString_(unit_accession, attributes, s_unit_accession);
      optionalAttributeAsString_(cv_ref, attributes, s_cv_ref);
      handleCVParam_(parent_parent_tag, parent_tag, attributeAsString_(attributes, s_accession), attributeAsString_(attributes, s_name), value, attributes, cv_ref, unit_accession);
		}

		if (tag_ == "MzIdentML")
		{
			// TODO handle version
			return;
		}

		if (tag_ == "Peptide")
		{
			// start new peptide
			actual_peptide_ = AASequence();
			// <Peptide id="peptide_1_1" sequenceMass="1341.72522" sequenceLength="12">
      //   <Modification location="0" residues="R" monoisotopicMassDelta="127.063324">
      //     <cvParam accession="UNIMOD:29" name="SMA" cvRef="UNIMOD" />
      //   </Modification>
      //   <peptideSequence>DAGTISGLNVLR</peptideSequence>
    	// </Peptide>
		}

		if (tag_ == "Modification")
		{
			
		}

		if (tag_ == "SpectrumIdentificationList")
		{
			
			return;
		}

		if (tag_ == "SpectrumIdentificationResult")
		{

			return;
		}

		if (tag_ == "SpectrumIdentificationItem")
		{
			//  <SpectrumIdentificationItem id="SII_1_1"  calculatedMassToCharge="670.86261" chargeState="2" experimentalMassToCharge="671.9" Peptide_ref="peptide_1_1" rank="1" passThreshold="true">
			// required attributes
			current_id_hit_.setId((attributeAsString_(attributes, "id")));
			current_id_hit_.setPassThreshold(asBool_(attributeAsString_(attributes, "passThreshold")));
			current_id_hit_.setRank(attributeAsInt_(attributes, "rank"));

			// optional attributes
			DoubleReal double_value(0);
			if (optionalAttributeAsDouble_(double_value, attributes, "calculatedMassToCharge"))
			{
				current_id_hit_.setCalculatedMassToCharge(double_value);
			}
			
			Int int_value(0);
			if (optionalAttributeAsInt_(int_value, attributes, "chargeState"))
			{
				current_id_hit_.setCharge(int_value);
			}

			if (optionalAttributeAsDouble_(double_value, attributes, "experimentalMassToCharge"))
			{
				current_id_hit_.setExperimentalMassToCharge(double_value);
			}

			if (optionalAttributeAsDouble_(double_value, attributes, "calculatedMassToCharge"))
			{
				current_id_hit_.setCalculatedMassToCharge(double_value);
			}
			
			String string_value("");
			if (optionalAttributeAsString_(string_value, attributes, "name"))
			{
				current_id_hit_.setName(string_value);
			}

			// TODO PeptideEvidence, pf:cvParam, pf:userParam, Fragmentation

			return;
		}

	}

	void MzIdentMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
	{
		if (tag_ == "Customizations")
		{
			String customizations = sm_.convert(chars);
			// TODO write customizations to Sofware
			return;
		}

		if (tag_ == "seq")
		{
			String seq = sm_.convert(chars);
			// TODO write to actual db sequence
			return;
		}

		if (tag_ == "peptideSequence")
		{
			String pep = sm_.convert(chars);
			//pep_sequences_[current] = current_sequence_;
			return;
		}

		// TODO any more?
	}

	void MzIdentMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		tag_ = sm_.convert(qname);
		open_tags_.pop_back();

		if (tag_ == "DataCollection")
		{
			return;
		}

		if (tag_ == "AnalysisData")
		{
			return;
		}

		if (tag_ == "ProteinDetectionList")
		{
			return;
		}

		if (tag_ == "SpectrumIdentificationList")
		{
			return;
		}

		if (tag_ == "SpectrumIdentificationResult")
		{
			return;
		}

		if (tag_ == "SpectrumIdentificationItem")
		{
			current_spectrum_id_.addHit(current_id_hit_);
			current_id_hit_ = IdentificationHit();
			return;
		}
	}

	void MzIdentMLHandler::handleCVParam_(const String& parent_parent_tag, const String& parent_tag, const String& accession, const String& name, const String& value, const xercesc::Attributes& attributes, const String& cv_ref, const String& unit_accession)	
	{
		if (parent_tag == "Modification")
		{
			if (cv_ref == "UNIMOD")
			{
				 //void ModificationsDB::searchModifications(set<const ResidueModification*>& mods, const String& origin, const String& name, ResidueModification::Term_Specificity term_spec) const
				set<const ResidueModification*> mods;
				Int loc = numeric_limits<Size>::max();
				if (optionalAttributeAsInt_(loc, attributes, "location"))
				{
					String uni_mod_id = accession.suffix(':');
					// TODO handle ambiguous residues
					String residues;
					if (optionalAttributeAsString_(residues, attributes, "residues"))
					{
						
					}
					if (loc == 0)
					{
        		ModificationsDB::getInstance()->searchTerminalModifications(mods, uni_mod_id, ResidueModification::N_TERM);
					}
					else if (loc == actual_peptide_.size())
					{
        		ModificationsDB::getInstance()->searchTerminalModifications(mods, uni_mod_id, ResidueModification::C_TERM);
					}
					else
					{
        		ModificationsDB::getInstance()->searchModifications(mods, residues, uni_mod_id, ResidueModification::ANYWHERE);
					}
				}
				else
				{
					warning(LOAD, "location of modification not defined!");
				}
			}
		}
	}

	void MzIdentMLHandler::writeTo(std::ostream& /*os*/)
	{
		// TODO
	}


	} //namespace Internal
} // namespace OpenMS



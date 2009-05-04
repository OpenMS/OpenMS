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
  	cv_.loadFromOBO("PI",File::find("/CV/psi-pi.obo"));
  }

  MzIdentMLHandler::MzIdentMLHandler(Identification& id, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			id_(&id),
			cid_(0)
  {
  	cv_.loadFromOBO("PI",File::find("/CV/psi-pi.obo"));
  }	

	MzIdentMLHandler::~MzIdentMLHandler()
	{
	}

	void MzIdentMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		tag_ = sm_.convert(qname);
		if (tag_ == "MzIdentML")
		{
			// TODO handle version
			return;
		}

		if (tag_ == "Peptide")
		{
			// start new peptide

			// <Peptide id="peptide_1_1" sequenceMass="1341.72522" sequenceLength="12">
      //   <Modification location="0" residues="R" monoisotopicMassDelta="127.063324">
      //     <pf:cvParam accession="UNIMOD:29" name="SMA" cvRef="UNIMOD" />
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
	
	void MzIdentMLHandler::writeTo(std::ostream& /*os*/)
	{
		// TODO
	}


	} //namespace Internal
} // namespace OpenMS



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

#include <OpenMS/FORMAT/HANDLERS/TraMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{

  TraMLHandler::TraMLHandler(const MRMExperiment& exp, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			exp_(0),
			cexp_(&exp)
  {
  	cv_.loadFromOBO("PI",File::find("/CV/psi-ms.obo"));
  }

  TraMLHandler::TraMLHandler(MRMExperiment& exp, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			exp_(&exp),
			cexp_(0)
  {
  	cv_.loadFromOBO("PI",File::find("/CV/psi-ms.obo"));
  }	

	TraMLHandler::~TraMLHandler()
	{
	}

	void TraMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		tag_ = sm_.convert(qname);
    open_tags_.push_back(tag_);

    static set<String> tags_to_ignore;
    if (tags_to_ignore.size() == 0)
    {
			tags_to_ignore.insert("TraML");
      tags_to_ignore.insert("contactList");
      tags_to_ignore.insert("cvList");
      tags_to_ignore.insert("instrumentList");
      tags_to_ignore.insert("softwareList");
      tags_to_ignore.insert("publicationList");
			tags_to_ignore.insert("proteinList");
    }

    // skip tags where nothing is to do
    if (tags_to_ignore.find(tag_) != tags_to_ignore.end())
    {
      return;
    }


    //determine parent tag
    String parent_tag;
    if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
    String parent_parent_tag;
    if (open_tags_.size()>2) parent_parent_tag = *(open_tags_.end()-3);

		if (tag_ == "cvParam")
		{
			String value = "";
      optionalAttributeAsString_(value, attributes, "value");
      String unit_accession = "";
      optionalAttributeAsString_(unit_accession, attributes, "unit_accession");
      handleCVParam_(parent_parent_tag, parent_tag, attributeAsString_(attributes, "accession"), attributeAsString_(attributes, "name"), value, unit_accession);
			return;
		}

		if (tag_ == "cv")
		{
			exp_->addCV(MRMExperiment::CV(attributeAsString_(attributes, "id"), attributeAsString_(attributes, "fullName"), attributeAsString_(attributes, "version"), attributeAsString_(attributes, "URI")));
			return;
		}

		if (tag_ == "contact")
		{
			MetaInfoInterface meta;
			String id = attributeAsString_(attributes, "id");
			meta.setMetaValue("id", id);
			actual_contact_ = meta;
			return;
		}

    if (tag_ == "publication")
    {
      MetaInfoInterface meta;
			String id = attributeAsString_(attributes, "id");
			meta.setMetaValue("id", id);
			actual_publication_ = meta;
      return;
    }

    if (tag_ == "instrument")
    {
     	MetaInfoInterface meta;
      String id = attributeAsString_(attributes, "id");
      meta.setMetaValue("id", id);
      actual_instrument_ = meta;
      return;
    }

    if (tag_ == "software")
    {
      actual_software_ = Software();
			actual_software_.setName(attributeAsString_(attributes, "id"));
			actual_software_.setVersion(attributeAsString_(attributes, "version"));
      return;
    }

    if (tag_ == "protein")
    {
      actual_protein_ = MRMExperiment::Protein();
			//<protein id="Q12149" accession="Q12149" name="Q12149" description="Heat shock protein 1, alpha" comment="">
			actual_protein_.id = attributeAsString_(attributes, "id");
			actual_protein_.accession = attributeAsString_(attributes, "accession");
			actual_protein_.name = attributeAsString_(attributes, "name");
			actual_protein_.description = attributeAsString_(attributes, "description");
			actual_protein_.comment = attributeAsString_(attributes, "comment");
      return;
    }

		if (tag_ == "peptide")
		{
			actual_peptide_ = MRMExperiment::Peptide();
			actual_peptide_.id = attributeAsString_(attributes, "id");
			actual_peptide_.protein_ref = attributeAsString_(attributes, "proteinRef");
			actual_peptide_.unmodified_sequence = attributeAsString_(attributes, "unmodifiedSequence");
			actual_peptide_.modified_sequence = attributeAsString_(attributes, "modifiedSequence");
			actual_peptide_.labeling_category = attributeAsString_(attributes, "labelingCategory");
			actual_peptide_.group_label = attributeAsString_(attributes, "groupLabel");
		}

		if (tag_ == "compound")
		{
			actual_compound_ = MRMExperiment::Compound();
			actual_compound_.id = attributeAsString_(attributes, "id");
		}

		if (tag_ == "retentionTime")
		{
			actual_rt_ = MRMExperiment::RetentionTime();
			DoubleReal local_retention_time(0);
			if (optionalAttributeAsDouble_(local_retention_time, attributes, "localRetentionTime"))
			{
				actual_rt_.local_retention_time = local_retention_time;
			}
			DoubleReal normalized_retention_time(0);
			if (optionalAttributeAsDouble_(normalized_retention_time, attributes, "normalizedRetentionTime"))
			{
				actual_rt_.normalized_retention_time = normalized_retention_time;
			}
			String normalization_standard;
			if (optionalAttributeAsString_(normalization_standard, attributes, "normalizationStandard"))
			{
				actual_rt_.normalization_standard = normalization_standard;
			}
			
			DoubleReal predicted_retention_time(0);
			if (optionalAttributeAsDouble_(predicted_retention_time, attributes, "predictedRetentionTime"))
			{
				actual_rt_.predicted_retention_time = predicted_retention_time;
			}

			String predicted_retention_time_software_ref;
			if (optionalAttributeAsString_(predicted_retention_time_software_ref, attributes, "predictedRetentionTimeSoftwareRef"))
			{
				actual_rt_.predicted_retention_time_software_ref = predicted_retention_time_software_ref;
			}
		}

		if (tag_ == "transition")
		{
			actual_transition_ = ReactionMonitoringTransition();
			String name;
			if (optionalAttributeAsString_(name, attributes, "name"))
			{
				actual_transition_.setName(name);
			}
			String peptide_ref;
			if (optionalAttributeAsString_(peptide_ref, attributes, "peptideRef"))
			{
				actual_transition_.setPeptideRef(peptide_ref);
			}
			String compound_ref;
			if (optionalAttributeAsString_(compound_ref, attributes, "compoundRef"))
			{
				actual_transition_.setCompoundRef(compound_ref);
			}
		}

		if (tag_ == "precursor")
		{
			actual_transition_.setPrecursorMZ(attributeAsDouble_(attributes, "mz"));
			Int charge(0);
			if (optionalAttributeAsInt_(charge, attributes, "charge"))
			{
				actual_transition_.setPrecursorCharge(charge);
			}
		}

    if (tag_ == "product")
    {
      actual_transition_.setProductMZ(attributeAsDouble_(attributes, "mz"));
      Int charge(0);
      if (optionalAttributeAsInt_(charge, attributes, "charge"))
      {
        actual_transition_.setProductCharge(charge);
      }
    }

		if (tag_ == "interpretation")
		{
			actual_interpretation_ = TransitionInterpretation();
			String product_series;
			if (optionalAttributeAsString_(product_series, attributes, "productSeries"))
			{
				actual_interpretation_.setProductSeries(product_series);
			}
			Int product_ordinal;
			if (optionalAttributeAsInt_(product_ordinal, attributes, "productOrdinal"))
			{
				actual_interpretation_.setProductOrdinal((Size)product_ordinal);
			}
			String product_adjustment;
			if (optionalAttributeAsString_(product_adjustment, attributes, "productAdjustment"))
			{
				actual_interpretation_.setProductAdjustment(product_adjustment);
			}
			DoubleReal mz_delta(0);
			if (optionalAttributeAsDouble_(mz_delta, attributes, "mzDelta"))
			{
				actual_interpretation_.setMZDelta(mz_delta);
			}
			String primary("false");
			if (optionalAttributeAsString_(primary, attributes, "primary"))
			{
				actual_interpretation_.setPrimary(DataValue(primary).toBool());
			}
		}


		cerr << "TraMLHandler: unknown tag opening: '" << tag_ << "'" << endl;
		return;
	}

	void TraMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
	{
		if (open_tags_.back() == "sequence")
		{
			String protein_sequence = sm_.convert(chars);
			actual_protein_.sequence = protein_sequence;
			return;
		}
		return;
	}

	void TraMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		tag_ = sm_.convert(qname);

		//determine parent tag
    String parent_tag;
    if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
    String parent_parent_tag;
    if (open_tags_.size()>2) parent_parent_tag = *(open_tags_.end()-3);

		open_tags_.pop_back();

		static set<String> tags_to_ignore;
		if (tags_to_ignore.size() == 0)
		{
			tags_to_ignore.insert("contactList");
			tags_to_ignore.insert("cvList");
			tags_to_ignore.insert("instrumentList");
			tags_to_ignore.insert("softwareList");
			tags_to_ignore.insert("publicationList");
			tags_to_ignore.insert("proteinList");
			tags_to_ignore.insert("precursor");
			tags_to_ignore.insert("product");
		}

		// skip tags where nothing is to do
		if (tags_to_ignore.find(tag_) != tags_to_ignore.end())
		{
			return;
		}

		if (tag_ == "cvParam")
		{
			return;
		}

		if (tag_ == "contact")
		{
			exp_->addContact(actual_contact_);
			return;
		}

		if (tag_ == "cv")
		{
			return;
		}

		if (tag_ == "instrument")
		{
			exp_->addInstrument(actual_instrument_);
			return;
		}

		if (tag_ == "publication")
		{
			exp_->addPublication(actual_publication_);
			return;
		}

		if (tag_ == "software")
		{
			exp_->addSoftware(actual_software_);
			return;
		}

		if (tag_ == "protein")
		{
			exp_->addProtein(actual_protein_);
			return;
		}

		if (tag_ == "retentionTime")
		{
			if (parent_tag == "peptide")
			{
				actual_peptide_.rts.push_back(actual_rt_);
				return;
			}
			if (parent_tag == "compound")
			{
				actual_compound_.rts.push_back(actual_rt_);
				return;
			}
			cerr << "TraMLHandler: tag 'retentionTime' no allowed at parent tag '" << parent_tag << "'" << endl;
		}

		if (tag_ == "peptide")
		{
			exp_->addPeptide(actual_peptide_);
		}

		if (tag_ == "compound")
		{
			exp_->addCompound(actual_compound_);
		}

		if (tag_ == "transition")
		{
			exp_->addTransition(actual_transition_);
		}
	
		if (tag_ == "interpretation")
		{
			actual_transition_.addInterpretation(actual_interpretation_);
		}
	
		cerr << "TraMLHandler: unknown tag closing: '" << tag_ << "'" << endl;
		return;
	}
	
  void TraMLHandler::writeTo(std::ostream& os)
  {
    const MRMExperiment& exp = *(cexp_);
    //logger_.startProgress(0,exp.size(),"storing mzML file");

    os  << "<?xml version=\"1.0\" encoding=\"UTF8\"?>" << endl;
    os  << "<TraML version=\"0.20\" xmlns=\"http://psi.hupo.org/ms/traml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/traml TraML0.2.xsd\">" << endl;
    //--------------------------------------------------------------------------------------------
    // CV list
    //--------------------------------------------------------------------------------------------
    os  << "  <cvList>" << endl;

		if (exp.getCVs().size() == 0)
		{
      os  << "    <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"unknown\" URI=\"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\"/>" << endl
        	<< "    <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"unknown\" URI=\"http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo\"/>" << endl;
		}
		else
		{
			for (vector<MRMExperiment::CV>::const_iterator it = exp.getCVs().begin(); it != exp.getCVs().end(); ++it)
			{
				os << "    <cv id=\"" << it->id << "\" fullName=\"" << it->fullname << "\" version=\"" << it->version << "\" URI=\"" << it->URI << "\"/>" << endl;
			}
		}
    os  << "  </cvList>" << endl;

		// contact list
		if (exp.getContacts().size() > 0)
		{
			os << "  <contactList>" << endl;
			for (vector<MetaInfoInterface>::const_iterator it = exp.getContacts().begin(); it != exp.getContacts().end(); ++it)
      {
        os << "    <contact id=\"" << (String)it->getMetaValue("id") << "\">" << endl;
        MetaInfoInterface meta = *it;
        meta.removeMetaValue("id");
        writeUserParam_(os, meta, 3);
        os << "    </contact>" << endl;
      }
      os << "  </contactList>" << endl;
		}

    // publication list
		if (exp.getPublications().size() > 0)
		{
			os << "  <publicationList>"  << endl;
			for (vector<MetaInfoInterface>::const_iterator it = exp.getPublications().begin(); it != exp.getPublications().end(); ++it)
			{	
				os << "    <publication id=\"" << (String)it->getMetaValue("id") << "\">" << endl;
				MetaInfoInterface meta = *it;
				meta.removeMetaValue("id");
				writeUserParam_(os, meta, 3);
				os << "    </publication>" << endl;
			}
			os << "  </publicationList>" << endl;
		}

    // instrument list
		if (exp.getInstruments().size() > 0)
		{
			os << "  <instrumentList>" << endl;
			for (vector<MetaInfoInterface>::const_iterator it = exp.getInstruments().begin(); it != exp.getInstruments().end(); ++it)
			{
				os << "    <instrument id=\"" << (String)it->getMetaValue("id") << "\">" << endl;
				MetaInfoInterface meta = *it;
				meta.removeMetaValue("id");
				writeUserParam_(os, meta, 3);
				os << "    </instrument>" << endl;
			}
			os << "  </instrumentList>" << endl;
		}

    // software list
		if (exp.getSoftware().size() > 0)
		{
			os << "  <softwareList>" << endl;
			for (vector<Software>::const_iterator it = exp.getSoftware().begin(); it != exp.getSoftware().end(); ++it)
			{
				os << "    <software id=\"" << it->getName() << "\" version=\"" << it->getVersion() << "\">" << endl;
				writeUserParam_(os, *it, 3);
				os << "    </software>" << endl;
			}
			os << "  </softwareList>" << endl;
		}

    // protein list
		if (exp.getProteins().size() > 0)
		{
			os << "  <proteinList>" << endl;
			for (vector<MRMExperiment::Protein>::const_iterator it = exp.getProteins().begin(); it != exp.getProteins().end(); ++it)
			{
				os << "    <protein id=\"" << it->id << "\" accession=\"" << it->accession << "\" name=\"" << it->name << "\" description=\"" << it->description << "\" comment=\"" << it->comment << "\" >" << endl;
				os << "      <sequence>" << it->sequence << "</sequence>" << endl;
				os << "    </protein>" << endl;
			}
			os << "  </proteinList>" << endl;
		}

    // compound list
		if (exp.getCompounds().size()  + exp.getPeptides().size() > 0)
		{
			os << "  <compoundList>" << endl;
			for (vector<MRMExperiment::Peptide>::const_iterator it = exp.getPeptides().begin(); it != exp.getPeptides().end(); ++it)
			{
				os << "    <peptide id=\"" << it->id << "\" proteinRef=\"" << it->protein_ref << "\" unmodifiedSequence=\"" << it->unmodified_sequence << "\" modifiedSequence=\"" << it->modified_sequence << "\" labelingCategory=\"" << it->labeling_category << "\" groupLabel=\"" << it->group_label << "\">" << endl;
				writeUserParams_(os, it->cvs, 3);
				
				for (vector<MRMExperiment::RetentionTime>::const_iterator rit = it->rts.begin(); rit != it->rts.end(); ++rit)
				{
					os << "       <retentionTime";
					if (rit->local_retention_time > 0)
					{
						os << " localRetentionTime=\"" << rit->local_retention_time << "\"";
					}
					if (rit->normalized_retention_time > 0)
					{
						os << " normalizedRetentionTime=\"" << rit->normalized_retention_time << "\"";
					}
					if (rit->normalization_standard != "")
					{
						os << " normalizationStandard=\"" << rit->normalization_standard << "\"";
					}
					if (rit->predicted_retention_time > 0)
					{
						os << " predictedRetentionTime=\"" << rit->predicted_retention_time << "\"";
					}
					if (rit->predicted_retention_time_software_ref != "")
					{
						os << " predictedRetentionTimeSoftwareRef=\"" << rit->predicted_retention_time_software_ref << "\"";
					}
					os << " >" << endl;
					writeUserParams_(os, rit->cvs, 5);
					os << "       </retentionTime>" << endl;
				}

				os << "    </peptide>" << endl;
			}

      for (vector<MRMExperiment::Compound>::const_iterator it = exp.getCompounds().begin(); it != exp.getCompounds().end(); ++it)
      {
        os << "    <compound id=\"" << it->id << "\">" << endl;
        writeUserParams_(os, it->cvs, 3);

				for (vector<MRMExperiment::RetentionTime>::const_iterator rit = it->rts.begin(); rit != it->rts.end(); ++rit)
        {
          os << "       <retentionTime";
          if (rit->local_retention_time > 0)
          {
            os << " localRetentionTime=\"" << rit->local_retention_time << "\"";
          }
          if (rit->normalized_retention_time > 0)
          {
            os << " normalizedRetentionTime=\"" << rit->normalized_retention_time << "\"";
          }
          if (rit->normalization_standard != "")
          {
            os << " normalizationStandard=\"" << rit->normalization_standard << "\"";
          }
          if (rit->predicted_retention_time > 0)
          {
            os << " predictedRetentionTime=\"" << rit->predicted_retention_time << "\"";
          }
          if (rit->predicted_retention_time_software_ref != "")
          {
            os << " predictedRetentionTimeSoftwareRef=\"" << rit->predicted_retention_time_software_ref << "\"";
          }
          os << " >" << endl;
          writeUserParams_(os, rit->cvs, 5);
          os << "       </retentionTime>" << endl;
        }

				os << "    </compound>" << endl;
			}

			os << "  </compoundList>" << endl;
		}

    // transition list
		if (exp.getTransitions().size() > 0)
		{
			/*
			os << "  <transitionList>" << endl;
			for (vector<ReactionMonitoringTransition>::const_iterator it = exp.getTransitions().begin(); it != exp.getTransitions().end(); ++it)
			{
				os << "    <transition";
				if (it->getName() != "")
				{
					os << " name=\"" << it->getName() << "\"";
				}

				if (it->getPeptideRef() != "")
				{
					os << " peptideRef=\"" << it->getPeptideRef() << "\"";
				}

				if (it->getCompoundRef() != "")
				{
					os << " compoundRef=\"" << it->getCompoundRef() << "\"";
				}
				os << " >" << endl;

				os << "      <precursor mz=\"" << it->getPrecursorMZ() << "\"";
				if (it->getPrecursorCharge() != numeric_limits<Int>::max())
				{
					os << " charge=\"" << it->getPrecursorCharge() << "\"";
				}
				os << " />" << endl;
			
				os << "      <product mz=\"" << it->getProductMZ() << "\"";
				if (it->getProductCharge() != numeric_limits<Int>::max())
				{
				 	os << "  charge=\"" << it->getProductCharge() << "\"";
				}
				os << " />" << endl;
				os << "      <interpretationList>" << endl;
				for (vector<TransitionInterpretation>::const_iterator tit = it->getInterpretations().begin(); tit != it->getInterpretations().end(); ++tit)
				{
					os << "        <interpretation";
					//<interpretation productSeries="y" productOrdinal="8" productAdjustment="" mzDelta="0.03" primary="true"/>
					if (tit->getProductSeries() != "")
					{
						os << " productSeries=\"" << tit->getProductSeries() << "\"";
					}
					if (tit->getProductOrdinal() != numeric_limits<Size>::max())
					{
						os << " productOrdinal=\"" << tit->getProductOrdinal() << "\"";
					}
					if (tit->getProductAdjustment() != "")
					{
						os << " productAdjustment=\"" << tit->getProductAdjustment() << "\"";
					}
					if (tit->getMZDelta() != numeric_limits<DoubleReal>::max())
					{
						os << " mzDelta=\"" << tit->getMZDelta() << "\"";
					}
					if (tit->getPrimary() == 1)
					{
						os << " primary=\"true\"";
					}

					if (tit->getPrimary() == 0)
					{
						os << " primary=\"false\"";
					}
					os << " >" << endl; 

					writeUserParams_(os, tit->getCVs(), 4);

					os << "        </interpretation>" << endl;
				}
				os << "      </interpretationList>" << endl;

				// TODO
				//os << "      <predictions>

				
				if (it->getConfigurations().size() > 0)
				{
					os << "      <configurationList>" << endl;
					for (vector<ReactionMonitoring::Configuration>::const_iterator cit = it->getConfigurations().begin(); cit != it->getConfigurations().end(); ++cit)
					{
						os << "       <configuration instrumentRef=\"" << cit->instrument_ref;
						if (cit->contact_ref != "")
						{
							os << " contactRef=\"" << cit->contact_ref;
						}
						os << " >" << endl;

						writeUserParams_(os, cit->cvs, 4);

						//for ()


						os << "       </configuration>" << endl;
					}
					os << "      </configurationList>" << endl;
				}
				

				os << "    </transition>" << endl;
			}
			os << "  </transitionList>" << endl;
			*/
		}


    os << "</TraML>" << endl;
    return;
  }

	void TraMLHandler::writeUserParam_(ostream& os, const MetaInfoInterface& meta, UInt indent) const
	{
		vector<String> keys;
		meta.getKeys(keys);
		for (vector<String>::const_iterator it = keys.begin(); it != keys.end(); ++it)
		{
			String key = *it;
			if (cv_.exists(key))
			{
				ControlledVocabulary::CVTerm term = cv_.getTerm(key);
				os << String(2 * indent, ' ') << "<cvParam cvRef=\"" << key.prefix(':') << "\" accession=\"" << key << "\" name=\"" << term.name << "\"";
				if (term.xref_type != ControlledVocabulary::CVTerm::NONE)
				{
					os << " value=\"" << (String)meta.getMetaValue(key) << "\"";
				}
				os << " />" << endl;
			}
			else
			{
				cerr << "TraMLHandler: unknown CV term '" << key << "' ignoring!" << endl;
			}
		}
	}

	void TraMLHandler::writeUserParams_(ostream& os, const vector<MetaInfoInterface>& meta, UInt indent) const
	{
		for (vector<MetaInfoInterface>::const_iterator it = meta.begin(); it != meta.end(); ++it)
		{
			writeUserParam_(os, *it, indent);
		}
	}


	void TraMLHandler::handleCVParam_(const String& /*parent_parent_tag*/, const String& parent_tag, const String& accession, const String& name, const String& value, const String& /*unit_accession*/)
	{
      //Error checks of CV values
      if (cv_.exists(accession))
      {
        const ControlledVocabulary::CVTerm& term = cv_.getTerm(accession);
        //obsolete CV terms
        if (term.obsolete)
        {
          warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
        }
        //check if term name and parsed name match
        String parsed_name = name;
        parsed_name.trim();
        String correct_name = term.name;
        correct_name.trim();
        if (parsed_name!=correct_name)
        {
          warning(LOAD, String("Name of CV term not correct: '") + term.id + " - " + parsed_name + "' should be '" + correct_name + "'");
        }
        if (term.obsolete)
        {
          warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
        //values used in wrong places and wrong value types
        if (value!="")
        {
          if (term.xref_type==ControlledVocabulary::CVTerm::NONE)
          {
            //Quality CV does not state value type :(
            if (!accession.hasPrefix("PATO:"))
            {
              warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must not have a value. The value is '" + value + "'.");
            }
          }
          else
          {
            switch(term.xref_type)
            {
              //string value can be anything
              case ControlledVocabulary::CVTerm::XSD_STRING:
                break;
              //int value => try casting
              case ControlledVocabulary::CVTerm::XSD_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NEGATIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_POSITIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NON_NEGATIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NON_POSITIVE_INTEGER:
                try
                {
                  value.toInt();
                }
                catch(Exception::ConversionError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must have an integer value. The value is '" + value + "'.");
                  return;
                }
                break;
              //double value => try casting
              case ControlledVocabulary::CVTerm::XSD_DECIMAL:
                try
                {
                  value.toDouble();
                }
                catch(Exception::ConversionError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must have a floating-point value. The value is '" + value + "'.");
                  return;
                }
                break;
              //date string => try conversion
              case ControlledVocabulary::CVTerm::XSD_DATE:
                try
                {
                  DateTime tmp;
                  tmp.set(value);
                }
                catch(Exception::ParseError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must be a valid date. The value is '" + value + "'.");
                  return;
                }
                break;
              default:
                warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' has the unknown value type '" + ControlledVocabulary::CVTerm::getXRefTypeName(term.xref_type) + "'.");
                break;
            }
          }
        }
        //no value, although there should be a numerical value
        else if (term.xref_type!=ControlledVocabulary::CVTerm::NONE && term.xref_type!=ControlledVocabulary::CVTerm::XSD_STRING)
        {
          warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' should have a numerical value. The value is '" + value + "'.");
          return;
        }
      }
		}


		// now handle the CVTerm and add it to the object
		if (parent_tag == "software")
		{
			actual_software_.setMetaValue(accession, value);
		}

		if (parent_tag == "publication")
		{
			actual_publication_.setMetaValue(accession, value);
		}

		if (parent_tag == "instrument")
		{
			actual_instrument_.setMetaValue(accession, value);
		}

		if (parent_tag == "contact")
		{
			actual_contact_.setMetaValue(accession, value);
		}

		if (parent_tag == "retentionTime")
		{
			MetaInfoInterface meta;
			meta.setMetaValue(accession, value);
			actual_rt_.cvs.push_back(meta);
		}
		
		if (parent_tag == "evidence")
		{
			MetaInfoInterface meta;
			meta.setMetaValue(accession, value);
			actual_peptide_.evidence.push_back(meta);
		}

		if (parent_tag == "peptide")
		{
			MetaInfoInterface meta;
			meta.setMetaValue(accession, value);
			actual_peptide_.cvs.push_back(meta);
		}

    if (parent_tag == "compound")
    {
      MetaInfoInterface meta;
      meta.setMetaValue(accession, value);
      actual_compound_.cvs.push_back(meta);
    }
	}

	} //namespace Internal
} // namespace OpenMS



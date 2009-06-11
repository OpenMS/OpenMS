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

		if (tag_ == "TraML")
		{
			// TODO handle version
			return;
		}

		if (tag_ == "cvList")
		{
			return;
		}

		if (tag_ == "cv")
		{
			exp_->addCV(MRMExperiment::CV(attributeAsString_(attributes, "id"), attributeAsString_(attributes, "fullName"), attributeAsString_(attributes, "version"), attributeAsString_(attributes, "URI")));
			return;
		}

		if (tag_ == "contactList")
		{
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

    if (tag_ == "publicationList")
    {
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

    if (tag_ == "instumentList")
    {
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

    if (tag_ == "softwareList")
    {
      // TODO
      return;
    }

    if (tag_ == "software")
    {
      // TODO
      return;
    }

    if (tag_ == "proteinList")
    {
      // TODO
      return;
    }

    if (tag_ == "protein")
    {
      // TODO
      return;
    }

		cerr << "TraMLHandler: unknown tag opening: '" << tag_ << "'" << endl;
		return;
	}

	void TraMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
	{
		if (open_tags_.back() == "protein")
		{
			String protein_sequence = sm_.convert(chars);
			return;
		}
		return;
	}

	void TraMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		tag_ = sm_.convert(qname);
		open_tags_.pop_back();

		if (tag_ == "cvParam")
		{
			return;
		}
		if (tag_ == "contactList")
		{
			return;
		}

		if (tag_ == "contact")
		{
			// TODO
			return;
		}

		if (tag_ == "cvList")
		{
			return;
		}
		
		if (tag_ == "cv")
		{
			return;
		}

		if (tag_ == "instrumentList")
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
		if (tag_ == "publicationList")
		{
			// nothing to do here
			return;
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

    // publication list
		if (exp.getPublications().size() > 0)
		{
			//os << "  <publicationList>"  << endl;
			//os << "  </publicationList>" << endl;
		}

    // instrument list
		if (exp.getInstruments().size() > 0)
		{
			//os << "  <instrumentList>" << endl;
			//os << "  </instrumentList>" << endl;
		}

    // software list
		if (exp.getSoftware().size() > 0)
		{
			//os << "  <softwareList>" << endl;
			//os << "  </softwareList>" << endl;
		}

    // protein list
		if (exp.getProteins().size() > 0)
		{
			os << "  <proteinList>" << endl;
			for (vector<MRMExperiment::Protein>::const_iterator it = exp.getProteins().begin(); it != exp.getProteins().end(); ++it)
			{
				os << "    <protein id=\"" << it->id << "\" accession=\"" << it->accession << "\" name=" << it->name << "\" description=\"" << it->description << "\" comment=\"" << it->comment << "\"" << endl;
				os << "      <sequence>" << endl;
				os << "        " << it->sequence << endl;
				os << "      </sequence>" << endl;
				os << "    </protein>" << endl;
			}
			os << "  </proteinList>" << endl;
		}

    // compound list
		if (exp.getCompounds().size()  + exp.getPeptides().size() > 0)
		{
			os << "  <compoundList>" << endl;
			for (vector<MRMExperiment::Compound>::const_iterator it = exp.getCompounds().begin(); it != exp.getCompounds().end(); ++it)
			{
				os << "    <compound id=\"" << it->id << "\">" << endl;
				writeUserParams_(os, it->cvs, 3);
        for (vector<MRMExperiment::RetentionTime>::const_iterator rit = it->rts.begin(); rit != it->rts.end(); ++rit)
        {
          os << "      <retentionTime localRetentionTime=\"" << rit->local_retention_time << "\" normalizedRetentionTime=\"" << rit->normalized_retention_time << "\" << normalizationStandard=\"" << rit->normalization_standard << "\" predictedRetentionTime=\"" << rit->predicted_retention_time << "\" predictedRetentionTimeSoftwareRef=\"" << rit->predicted_retention_time_software_ref << "\" >" << endl;
          writeUserParams_(os, rit->cvs, 5);
          os << "      </retentionTime>" << endl;
        }
				os << "    </compound>" << endl;
			}

			for (vector<MRMExperiment::Peptide>::const_iterator it = exp.getPeptides().begin(); it != exp.getPeptides().end(); ++it)
			{
				os << "    <peptide id=\"" << it->id << "\" proteinRef=\"" << it->protein_ref << "\" unmodifiedSequence=\"" << it->unmodified_sequence << "\" modifiedSequence=\"" << it->modified_sequence << "\" labelingCategory=\"" << it->labeling_category << "\" groupLabel=\"" << it->group_label << "\">" << endl;
				writeUserParams_(os, it->cvs, 3);
				
				for (vector<MRMExperiment::RetentionTime>::const_iterator rit = it->rts.begin(); rit != it->rts.end(); ++rit)
				{
					os << "      <retentionTime localRetentionTime=\"" << rit->local_retention_time << "\" normalizedRetentionTime=\"" << rit->normalized_retention_time << "\" << normalizationStandard=\"" << rit->normalization_standard << "\" predictedRetentionTime=\"" << rit->predicted_retention_time << "\" predictedRetentionTimeSoftwareRef=\"" << rit->predicted_retention_time_software_ref << "\" >" << endl;
					writeUserParams_(os, rit->cvs, 5);
					os << "      </retentionTime>" << endl;
				}

				os << "    </peptide>" << endl;
			}
			os << "  </compoundList>" << endl;
		}

    // transition list
		if (exp.getTransitions().size() > 0)
		{
			os << "  <transitionList>" << endl;
			// TODO
			os << "  </transitionList>" << endl;
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
				os << String('\t', indent) << "<cvParam cvRef=\"" << key.prefix(':') << "\" accession=\"" << key << "\" name=\"" << term.name;
				if (term.xref_type != ControlledVocabulary::CVTerm::NONE)
				{
					os << "\" value=\"" << (String)meta.getMetaValue(key) << "\"";
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


		
	}

	} //namespace Internal
} // namespace OpenMS



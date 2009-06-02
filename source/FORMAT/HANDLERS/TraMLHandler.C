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
		if (tag_ == "TraML")
		{
			// TODO handle version
			return;
		}

		if (tag_ == "cvList")
		{
			// TODO
			return;
		}

		if (tag_ == "cv")
		{
			// TODO
			return;
		}

		if (tag_ == "contactList")
		{
			// TODO
			return;
		}

		if (tag_ == "contact")
		{
			// TODO 
			return;
		}

    if (tag_ == "publicationList")
    {
      // TODO
      return;
    }

    if (tag_ == "publication")
    {
      // TODO
      return;
    }

    if (tag_ == "instumentList")
    {
      // TODO
      return;
    }

    if (tag_ == "instrument")
    {
      // TODO
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

		return;
	}

	void TraMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
	{
		if (tag_ == "protein")
		{
			String protein_sequence = sm_.convert(chars);
			return;
		}
		return;
	}

	void TraMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		tag_ = sm_.convert(qname);
		if (tag_ == "contactList")
		{
			return;
		}

		if (tag_ == "contact")
		{
			return;
		}

		if (tag_ == "ProteinDetectionList")
		{
			return;
		}

		if (tag_ == "SpectrumMRMExperimentList")
		{
			return;
		}

		if (tag_ == "SpectrumMRMExperimentResult")
		{
			return;
		}

		if (tag_ == "SpectrumMRMExperimentItem")
		{
			return;
		}
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
				os << "    <cv id=\"" << it->id << "\" fullname=\"" << it->fullname << "\"version=\"" << it->version << "\" URI=\"" << it->URI << "\"/>" << endl;
			}
		}
    os  << "  </cvList>" << endl;

    // publication list
		if (exp.getPublications().size() > 0)
		{
			os << "  <publicationList>"  << endl;
			os << "  </publicationList>" << endl;
		}

    // instrument list
		if (exp.getInstruments().size() > 0)
		{
			os << "  <instrumentList>" << endl;
			os << "  </instrumentList>" << endl;
		}

    // software list
		if (exp.getSoftware().size() > 0)
		{
			os << "  <softwareList>" << endl;
			os << "  </softwareList>" << endl;
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
		/*
		if (exp.getCompounts().size() > 0)
		{
			os << "  <compoundList>" << endl;
			os << "  </compoundList>" << endl;
		}
		*/

    // transition list
		if (exp.getTransitions().size() > 0)
		{
			os << "  <transitionList>" << endl;
			os << "  </transitionList>" << endl;
		}


    os << "</TraML>" << endl;
    return;
  }



	} //namespace Internal
} // namespace OpenMS



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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FileConverter FileConverter

	@brief Converts between different MS file formats.

	This converter tries to determine the file type from the file extension or from the first few lines
	of the file. If file type determination is not possible, you have to give the input or output file type explicitly.

	During some conversion operations information is lost, e.g. when converting featureXML to mzData.
	In these cases a warning is shown.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FileConverter.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileConverter
	: public TOPPBase
{
 public:
	TOPPFileConverter()
		: TOPPBase("FileConverter","Converts between different MS file formats.")
	{

	}

 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file ");
		setValidFormats_("in",StringList::create("mzData,mzXML,mzML,DTA,DTA2D,cdf,mgf,featureXML,consensusXML,ms2"));
		registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
		setValidStrings_("in_type",StringList::create("mzData,mzXML,mzML,DTA,DTA2D,cdf,mgf,featureXML,consensusXML,ms2"));

		registerOutputFile_("out","<file>","","output file ");
		setValidFormats_("out",StringList::create("mzData,mzXML,mzML,DTA2D,mgf,featureXML"));
		registerStringOption_("out_type", "<type>", "", "output file type -- default: determined from file extension or content\n", false);
		setValidStrings_("out_type",StringList::create("mzData,mzXML,mzML,DTA2D,mgf,featureXML"));
	}

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------

		//input file names
		String in = getStringOption_("in");

		//input file type
		FileHandler fh;
		FileTypes::Type in_type = fh.nameToType(getStringOption_("in_type"));

		if (in_type==FileTypes::UNKNOWN)
		{
			in_type = fh.getType(in);
			writeDebug_(String("Input file type: ") + fh.typeToName(in_type), 2);
		}

		if (in_type==FileTypes::UNKNOWN)
		{
			writeLog_("Error: Could not determine input file type!");
			return PARSE_ERROR;
		}


		//output file names and types
		String out = getStringOption_("out");
		FileTypes::Type out_type = fh.nameToType(getStringOption_("out_type"));

		if (out_type==FileTypes::UNKNOWN)
		{
			out_type = fh.getTypeByFileName(out);
		}

		if (out_type==FileTypes::UNKNOWN)
		{
			writeLog_("Error: Could not determine output file type!");
			return PARSE_ERROR;
		}

		writeDebug_(String("Output file type: ") + fh.typeToName(out_type), 1);

		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
		typedef MSExperiment< Peak1D > MSExperimentType;
		MSExperimentType exp;

		typedef MSExperimentType::SpectrumType SpectrumType;

		typedef FeatureMap<> FeatureMapType;

    FeatureMapType fm;
    ConsensusMap cm;

		vector<ProteinIdentification> prot_ids;
		vector<PeptideIdentification> pep_ids;

		writeDebug_(String("Loading input file"), 1);

		if (in_type == FileTypes::FEATUREXML)
		{
      FeatureXMLFile().load(in,fm);
      fm.sortByPosition();
		  if ( out_type != FileTypes::FEATUREXML )
		  {
	      // You will lose information and waste memory. Enough reasons to issue a warning!
		    writeLog_("Warning: Converting features to peaks. You will lose information!");
	      exp.set2DData(fm);
		  }
		}
		else if (in_type == FileTypes::CONSENSUSXML)
		{
			ConsensusXMLFile().load(in,cm);
			cm.sortByPosition();
      if ( out_type != FileTypes::FEATUREXML )
      {
        // You you will lose information and waste memory. Enough reasons to issue a warning!
        writeLog_("Warning: Converting consensus features to peaks. You will lose information!");
        exp.set2DData(cm);
      }
		}
		else
		{
			fh.loadExperiment(in,exp,in_type,log_type_);
		}

		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------

		writeDebug_(String("Writing output file"), 1);
		if (out_type == FileTypes::MZML)
		{
			//add data processing entry
			addDataProcessing_(exp, getProcessingInfo_(DataProcessing::CONVERSION_MZML));

			MzMLFile f;
			f.setLogType(log_type_);
			ChromatogramTools().convertSpectraToChromatograms(exp, true);
			f.store(out,exp);
		}
		else if (out_type == FileTypes::MZDATA)
		{
			//annotate output with data processing info
			addDataProcessing_(exp, getProcessingInfo_(DataProcessing::CONVERSION_MZDATA));
			
			MzDataFile f;
			f.setLogType(log_type_);
			ChromatogramTools().convertChromatogramsToSpectra<MSExperimentType>(exp);
			f.store(out,exp);
		}
		else if (out_type == FileTypes::MZXML)
		{
			//annotate output with data processing info
			addDataProcessing_(exp, getProcessingInfo_(DataProcessing::CONVERSION_MZXML));
			
			MzXMLFile f;
			f.setLogType(log_type_);
			ChromatogramTools().convertChromatogramsToSpectra<MSExperimentType>(exp);
			f.store(out,exp);
		}
		else if (out_type == FileTypes::DTA2D)
		{
			//add data processing entry
			addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));
			
			DTA2DFile f;
			f.setLogType(log_type_);
			ChromatogramTools().convertChromatogramsToSpectra<MSExperimentType>(exp);
			f.store(out,exp);
		}
		else if (out_type == FileTypes::FEATUREXML)
		{
		  if ( in_type == FileTypes::FEATUREXML )
		  {
		    fm.applyMemberFunction(&UniqueIdInterface::setUniqueId);

		    //add data processing entry
        addDataProcessing_(fm, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

	      FeatureXMLFile().store(out,fm);
		  }
		  else if ( in_type == FileTypes::CONSENSUSXML )
		  {
		    fm.resize(cm.size());
		    // fm.MetaInfoInterface::operator=(cm); // not available ...
		    fm.DocumentIdentifier::operator=(cm);
		    fm.UniqueIdInterface::operator=(cm);
		    fm.setProteinIdentifications(cm.getProteinIdentifications());
        fm.setUnassignedPeptideIdentifications(cm.getUnassignedPeptideIdentifications());
		    for ( Size i = 0; i < cm.size(); ++i )
		    {
		      Feature & f = fm[i];
		      const ConsensusFeature & c = cm[i];
		      f.RichPeak2D::operator=(c);
		      f.setCharge(c.getCharge());
		      f.setOverallQuality(c.getQuality());
		      f.setPeptideIdentifications(c.getPeptideIdentifications());
		    }

		    // TODO Discuss this: Arguably we do NOT want to assign new unique ids?
		    // fm.applyMemberFunction(&UniqueIdInterface::setUniqueId);

        //add data processing entry
        addDataProcessing_(fm, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

		    FeatureXMLFile().store(out,fm);
		  }
		  else
		  {
        // The feature specific information is only defaulted. Enough reasons to issue a warning!
        writeLog_("Warning: Converting peaks to features results in incomplete features!");
        FeatureMapType feature_map;
        feature_map.reserve(exp.getSize());
        typedef FeatureMapType::FeatureType FeatureType;
        FeatureType feature;
        feature.setQuality(0,1); // override default
        feature.setQuality(1,1); // override default
        feature.setOverallQuality(1); // override default
        for ( MSExperimentType::ConstIterator spec_iter = exp.begin();
              spec_iter != exp.end();
              ++spec_iter
            )
        {
          feature.setRT(spec_iter->getRT());
          for ( SpectrumType::ConstIterator peak1_iter = spec_iter->begin();
                peak1_iter != spec_iter->end();
                ++peak1_iter
              )
          {
            feature.setMZ(peak1_iter->getMZ());
            feature.setIntensity(peak1_iter->getIntensity());
            feature.setUniqueId();
            feature_map.push_back(feature);
          }
        }
        feature_map.updateRanges();

        //add data processing entry
        addDataProcessing_(feature_map, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

        FeatureXMLFile().store(out,feature_map);
		  }
		}
		else if (out_type == FileTypes::MGF)
		{
			//add data processing entry
			addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));
				
			MascotGenericFile f;
			Param p(f.getParameters());
			p.setValue("peaklists_only", "true");
			f.setParameters(p);
			f.store(out, exp);
		}
		else
		{
			writeLog_("Unknown output file type given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}

		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TOPPFileConverter tool;
	return tool.main(argc,argv);
}

/// @endcond

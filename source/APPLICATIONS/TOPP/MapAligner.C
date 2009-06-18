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
// $Maintainer: Marc Sturm, Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>


#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_MapAligner MapAligner

	@brief Corrects retention time distortions between maps.

	This tool provides several different algorithms to correct for retention time shifts
	and distortions.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_MapAligner.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAligner
  : public TOPPBase
{

public:
	TOPPMapAligner()
		: TOPPBase("MapAligner","Corrects retention time distortions between maps.")
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		registerInputFileList_("in","<files>",StringList(),"input files separated by blanks",true);
		setValidFormats_("in",StringList::create("mzML,featureXML,idXML"));
		registerOutputFileList_("out","<files>",StringList(),"output files separated by blanks",false);
		setValidFormats_("out",StringList::create("mzML,featureXML,idXML"));
		registerOutputFileList_("transformations","<files>",StringList(),"transformation output files separated by blanks",false);
		registerStringOption_("type","<name>","","Map alignment algorithm type",true);
		setValidStrings_("type", getToolList()[toolName_()] );

		// TODO  Remove this hack when StringList when become available in INIFileEditor.
		registerInputFileList_("given_transformations","<files>",StringList(),"given transformations separated by blanks. [This is a workaround used by algorithm type apply_given_trafo until StringList is supported by INIFileEditor.]",false);
		setValidFormats_("given_transformations",StringList::create("trafoXML"));

    addEmptyLine_();
		addText_("This tool takes N input files, aligns them and writes them to the output files.");

		registerSubsection_("algorithm","Algorithm parameters section");
	}

	Param getSubsectionDefaults_(const String& /* section */ ) const
	{
		String type = getStringOption_("type");
		MapAlignmentAlgorithm* algo = Factory<MapAlignmentAlgorithm>::create(type);
		Param tmp = algo->getParameters();

		// TODO  Remove this hack when StringList when become available in INIFileEditor.
		if ( type == "apply_given_trafo")
		{
			tmp.setValue("transformations",getStringList_("given_transformations"));
		}

		delete algo;
		return tmp;
	}

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		StringList ins = getStringList_("in");

		StringList outs = getStringList_("out");

		StringList trafos = getStringList_("transformations");

		String type = getStringOption_("type");

		ProgressLogger progresslogger;
		progresslogger.setLogType(log_type_);

		//-------------------------------------------------------------
		// check for valid input
		//-------------------------------------------------------------
		//check whether the number of input files equals the number of output files
		if (ins.size()!=outs.size())
		{
			writeLog_("Error: The number of input and output files has to be equal!");
			return ILLEGAL_PARAMETERS;
		}
		//check whether the number of input files equals the number of output files
		if (trafos.size()!=0 && ins.size()!=trafos.size())
		{
			writeLog_("Error: The number of input and transformation files has to be equal!");
			return ILLEGAL_PARAMETERS;
		}
		//check whether all input files have the same type (this type is used to store the output type too)
		FileTypes::Type in_type = FileHandler::getType(ins[0]);
		for (Size i=1;i<ins.size();++i)
		{
			if (FileHandler::getType(ins[i])!=in_type)
			{
				writeLog_("Error: All input files have to be in the same format!");
				return ILLEGAL_PARAMETERS;
			}
		}

    //-------------------------------------------------------------
    // set up alignment algorithm
    //-------------------------------------------------------------
    MapAlignmentAlgorithm* alignment = Factory<MapAlignmentAlgorithm>::create(type);
		Param alignment_param = getParam_().copy("algorithm:", true);

		// TODO  Remove this hack when StringList when become available in INIFileEditor.
		if (type == "apply_given_trafo")
		{
			alignment_param.setValue("transformations",getStringList_("given_transformations"));
		}

		writeDebug_("Used alignment parameters", alignment_param, 3);
		alignment->setParameters(alignment_param);
		alignment->setLogType(log_type_);


    //-------------------------------------------------------------
    // perform peak alignment
    //-------------------------------------------------------------
		std::vector<TransformationDescription> transformations;
		if (in_type == FileTypes::MZML)
		{
			// load input
			std::vector< MSExperiment<> > peak_maps(ins.size());
			MzMLFile f;
			f.setLogType(log_type_);
			for (Size i=0; i<ins.size(); ++i)
			{
		    f.load(ins[i], peak_maps[i]);
			}

			// try to align
			try
			{
				alignment->alignPeakMaps(peak_maps,transformations);
			}
			catch (Exception::NotImplemented&)
			{
				writeLog_("Error: The algorithm '" + type + "' cannot be used for peak data!");
				return INTERNAL_ERROR;
			}

			// write output
			progresslogger.startProgress(0,outs.size(),"writing output files");
			for (Size i=0; i<outs.size(); ++i)
			{
				progresslogger.setProgress(i);
				
				//annotate output with data processing info
				addDataProcessing_(peak_maps[i], getProcessingInfo_(DataProcessing::ALIGNMENT));

		    f.store(outs[i], peak_maps[i]);
			}
			progresslogger.endProgress();
		}
    //-------------------------------------------------------------
    // perform feature alignment
    //-------------------------------------------------------------
		else if (in_type == FileTypes::FEATUREXML)
		{
			// load input
			std::vector< FeatureMap<> > feat_maps(ins.size());
			FeatureXMLFile f;
			// f.setLogType(log_type_); // TODO
			progresslogger.startProgress(0,ins.size(),"loading input files");
			for (Size i=0; i<ins.size(); ++i)
			{
				progresslogger.setProgress(i);
		    f.load(ins[i], feat_maps[i]);
			}
			progresslogger.endProgress();

			// try to align
			try
			{
				alignment->alignFeatureMaps(feat_maps,transformations);
			}
			catch (Exception::NotImplemented&)
			{
				writeLog_("Error: The algorithm '" + type + "' cannot be used for feature data!");
				return INTERNAL_ERROR;
			}

			// write output
			progresslogger.startProgress(0,outs.size(),"writing output files");
			for (Size i=0; i<outs.size(); ++i)
			{
				progresslogger.setProgress(i);
				
				//annotate output with data processing info
				addDataProcessing_(feat_maps[i], getProcessingInfo_(DataProcessing::ALIGNMENT));

		    f.store(outs[i], feat_maps[i]);
			}
			progresslogger.endProgress();
		}
    //-------------------------------------------------------------
    // perform peptide alignment
    //-------------------------------------------------------------
		else if (in_type == FileTypes::IDXML)
		{
			// load input
			std::vector<  std::vector<ProteinIdentification> > protein_ids_vec(ins.size());
			std::vector<  std::vector<PeptideIdentification> > peptide_ids_vec(ins.size());

			IdXMLFile f;
			// f.setLogType_(log_type_);

			progresslogger.startProgress(0,ins.size(),"loading input files");
			for (Size i=0; i<ins.size(); ++i)
			{
				progresslogger.setProgress(i);
				String document_id;
		    f.load( ins[i], protein_ids_vec[i], peptide_ids_vec[i], document_id);
			}
			progresslogger.endProgress();

			// try to align
			try
			{
				alignment->alignPeptideIdentifications(peptide_ids_vec,transformations);
			}
			catch (Exception::NotImplemented&)
			{
				writeLog_("Error: The algorithm '" + type + "' cannot be used for peptide data!");
				return INTERNAL_ERROR;
			}

			// write output
			progresslogger.startProgress(0,outs.size(),"writing output files");
			for (Size i=0; i<outs.size(); ++i)
			{
				progresslogger.setProgress(i);
		    f.store( outs[i], protein_ids_vec[i], peptide_ids_vec[i] );
			}
			progresslogger.endProgress();
		}
		else
		{
			// TODO can this really happen? I think it is tested above. Otherwise
			// throw an appropriate exception?
			return ILLEGAL_PARAMETERS;
		}

		delete alignment;

		if (trafos.size()!=0)
		{
			for (Size i=0; i<transformations.size(); ++i)
			{
				TransformationXMLFile().store(trafos[i],transformations[i]);
			}
		}

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
  TOPPMapAligner tool;
  return tool.main(argc,argv);
}

/// @endcond

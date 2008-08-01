// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>


#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page MapAligner MapAligner
 
	 @brief Corrects retention time distortions between maps.
	 
	 This tool provides several different algorithms to correct for retention time shifts
	 and distortions.
	 
	 @todo write map alignment algorithm that takes TrafoXML and applies it (Clemens)
	 
	 @ingroup TOPP
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
		registerStringOption_("in","<files>","","Comma-separated list of input file names in FeatureXML or mzData format",true);
		registerStringOption_("out","<files>","","Comma-separated list of output file names in FeatureXML or mzData format");
		registerStringOption_("transformations","<files>","","Comma-separated list of output files for transformations",false);
		registerStringOption_("type","<name>","","Map alignment algorithm type",true);
		setValidStrings_("type",Factory<MapAlignmentAlgorithm>::registeredProducts());
    
    addEmptyLine_();
		addText_("This tool takes N input files, aligns them and writes them to the output files.");
    
		registerSubsection_("algorithm","Algorithm parameters section");
	}
	
	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		String type = getStringOption_("type");
		MapAlignmentAlgorithm* alignment = Factory<MapAlignmentAlgorithm>::create(type);
		Param tmp = alignment->getParameters();
		delete alignment;
		return tmp;
	}   

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		StringList ins;
		String in = getStringOption_("in");
		in.split(',',ins);
		if (ins.size()==0) ins.push_back(in);

		StringList outs;
		String out = getStringOption_("out");
		out.split(',',outs);
		if (outs.size()==0) outs.push_back(out);

		StringList trafos;		
		if (setByUser_("transformations") && getStringOption_("transformations")!="")
		{
			String trafo = getStringOption_("transformations");
			trafo.split(',',trafos);
			if (trafos.size()==0) trafos.push_back(trafo);
		}

		String type = getStringOption_("type");
		
		//-------------------------------------------------------------
		// check for valid input
		//-------------------------------------------------------------
		//check if the numer of input files equals the number of output files
		if (ins.size()!=outs.size())
		{
			writeLog_("Error: The number of input and output files has to be equal!");
			return ILLEGAL_PARAMETERS;
		}
		//check if the numer of input files equals the number of output files
		if (trafos.size()!=0 && ins.size()!=trafos.size())
		{
			writeLog_("Error: The number of input and transformation files has to be equal!");
			return ILLEGAL_PARAMETERS;
		}
		//check if all input files have the same type (this type is used to store the output type too)		
		FileHandler::Type in_type = FileHandler::getType(ins[0]);
		for (UInt i=1;i<ins.size();++i)
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
		writeDebug_("Used alignment parameters", alignment_param, 3);
		alignment->setParameters(alignment_param);

    //-------------------------------------------------------------
    // perform peak alignment
    //-------------------------------------------------------------
		std::vector<TransformationDescription> transformations;
		if (in_type == FileHandler::MZDATA)
		{
			//load input
			std::vector< MSExperiment<> > peak_maps(ins.size());
			MzDataFile f;
			f.setLogType(log_type_);
			for (UInt i=0; i<ins.size(); ++i)
			{		 		
		    f.load(ins[i], peak_maps[i]);
			}
			
			//try to align
			try
			{
				alignment->alignPeakMaps(peak_maps,transformations);
			}
			catch (Exception::NotImplemented&)
			{
				writeLog_("Error: The algorithm '" + type + "' can only be used for feature data!");
				return INTERNAL_ERROR;
			}
			
			//write output
			for (UInt i=0; i<outs.size(); ++i)
			{		 		
		    f.store(outs[i], peak_maps[i]);
			}
		}
    //-------------------------------------------------------------
    // perform feature alignment
    //-------------------------------------------------------------
		else
		{
			//load input
			std::vector< FeatureMap<> > feat_maps(ins.size());
			FeatureXMLFile f;
			for (UInt i=0; i<ins.size(); ++i)
			{		 		
		    f.load(ins[i], feat_maps[i]);
			}

			//try to align
			try
			{
				alignment->alignFeatureMaps(feat_maps,transformations);
			}
			catch (Exception::NotImplemented&)
			{
				writeLog_("Error: The algorithm '" + type + "' can only be used for peak data!");
				return INTERNAL_ERROR;
			}
			
			//write output
			for (UInt i=0; i<outs.size(); ++i)
			{		 		
		    f.store(outs[i], feat_maps[i]);
			}
		}
		
		delete alignment;
		
		if (trafos.size()!=0)
		{
			for (UInt i=0; i<transformations.size(); ++i)
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

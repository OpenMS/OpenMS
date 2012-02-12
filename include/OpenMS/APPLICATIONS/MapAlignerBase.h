// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm, Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
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
// Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_MapAlignerBase MapAlignerBase

	@brief Base class for different MapAligner TOPP tools.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignerBase
  : public TOPPBase
{

public:
	TOPPMapAlignerBase(String name, String description, bool official=true)
		: TOPPBase(name, description, official)
	{
	}

	// "public" so it can be used in DefaultParamHandlerDocumenter to get docu
	static Param getModelDefaults(const String& default_model)
	{
		Param params;
		params.setValue("type", default_model, "Type of model");
		// TODO: avoid referring to each TransformationModel subclass explicitly
		StringList model_types = StringList::create("linear,b_spline,interpolated");
		if (!model_types.contains(default_model))
		{
			model_types.insert(model_types.begin(), default_model);
		}
		params.setValidStrings("type", model_types);

		Param model_params;
		TransformationModelLinear::getDefaultParameters(model_params);
		params.insert("linear:", model_params);
		params.setSectionDescription("linear", "Parameters for 'linear' model");
		TransformationModelBSpline::getDefaultParameters(model_params);
		params.insert("b_spline:", model_params);
		params.setSectionDescription("b_spline", "Parameters for 'b_spline' model");
		TransformationModelInterpolated::getDefaultParameters(model_params);
		// "polynomial" interpolation is not suitable for RT data, so remove it:
		const Param::ParamEntry& entry = 
			model_params.getEntry("interpolation_type");
		StringList interpolation_types = entry.valid_strings;
		StringList::Iterator pos = find(interpolation_types.begin(), 
																		interpolation_types.end(), "polynomial"); 
		interpolation_types.erase(pos); 
		model_params.setValidStrings("interpolation_type", interpolation_types);
		params.insert("interpolated:", model_params);
		params.setSectionDescription("interpolated", 
																 "Parameters for 'interpolated' model");
		return params;
	}

protected:
	void registerOptionsAndFlags_(const String& file_formats)
	{
		registerInputFileList_("in", "<files>", StringList(), "Input files separated by blanks (all must have the same file type)", true);
		setValidFormats_("in", StringList::create(file_formats));
		registerOutputFileList_("out", "<files>", StringList(), "Output files separated by blanks", false);
		setValidFormats_("out", StringList::create(file_formats));
		registerOutputFileList_("trafo_out", "<files>", StringList(), "Transformation output files separated by blanks", false);
		setValidFormats_("trafo_out", StringList::create("trafoXML"));
    addEmptyLine_();
		addText_("This tool takes a number of input files, aligns them and writes the results to the output files.");
		addText_("Either 'out' or 'trafo_out' has to be provided. They can be used together.");
	}

	void handleReference_(MapAlignmentAlgorithm* alignment)
	{
		// note: this function is in the base class to avoid code duplication, but
		// it only makes sense for some derived classes - don't call the function
		// in a class that doesn't support a reference!

		// check reference parameters:
		Size reference_index = getIntOption_("reference:index");
		String reference_file = getStringOption_("reference:file");
		if (reference_index > getStringList_("in").size())
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'reference:index' must not be higher than the number of input files");
		}
		if (reference_index && !reference_file.empty())
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'reference:index' and 'reference:file' cannot be used together");
		}

		// pass the reference parameters on to the algorithm:
		alignment->setReference(reference_index, reference_file);
	}

	ExitCodes commonMain_(MapAlignmentAlgorithm* alignment)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		StringList ins = getStringList_("in");
		StringList outs = getStringList_("out");
		StringList trafos = getStringList_("trafo_out");
		Param model_params = getParam_().copy("model:", true);
		String model_type = model_params.getValue("type");
		model_params = model_params.copy(model_type + ":", true);

		ProgressLogger progresslogger;
		progresslogger.setLogType(log_type_);

		//-------------------------------------------------------------
		// check for valid input
		//-------------------------------------------------------------
		// check whether some kind of output file is given:
		if (outs.empty() && trafos.empty())
		{
			writeLog_("Error: Either data output or transformation output files have to be provided!");
			return ILLEGAL_PARAMETERS;
		}
		// check whether number of input files equals number of output files:
		if (!outs.empty() && (ins.size() != outs.size()))
		{
			writeLog_("Error: The number of input and output files has to be equal!");
			return ILLEGAL_PARAMETERS;
		}
		if (!trafos.empty() && (ins.size() != trafos.size()))
		{
			writeLog_("Error: The number of input and transformation output files has to be equal!");
			return ILLEGAL_PARAMETERS;
		}
		// check whether all input files have the same type (this type is used to store the output type too):
		FileTypes::Type in_type = FileHandler::getType(ins[0]);
		for (Size i = 1; i < ins.size(); ++i)
		{
			if (FileHandler::getType(ins[i]) != in_type)
			{
				writeLog_("Error: All input files have to be in the same format!");
				return ILLEGAL_PARAMETERS;
			}
		}

    //-------------------------------------------------------------
    // set up alignment algorithm
    //-------------------------------------------------------------
		Param alignment_param = getParam_().copy("algorithm:", true);

		writeDebug_("Used alignment parameters", alignment_param, 3);
		alignment->setParameters(alignment_param);
		alignment->setLogType(log_type_);

		std::vector<TransformationDescription> transformations;

    //-------------------------------------------------------------
    // perform peak alignment
    //-------------------------------------------------------------
		if (in_type == FileTypes::MZML)
		{
			// load input
			std::vector< MSExperiment<> > peak_maps(ins.size());
			MzMLFile f;
			f.setLogType(log_type_);
			for (Size i = 0; i < ins.size(); ++i)
			{
		    f.load(ins[i], peak_maps[i]);
			}

			// try to align
			try
			{
				alignment->alignPeakMaps(peak_maps, transformations);
			}
			catch (Exception::NotImplemented&)
			{
				writeLog_("Error: The algorithm '" + alignment->getName() + "' cannot be used for peak data!");
				return INTERNAL_ERROR;
			}
			if (model_type != "none")
			{
				alignment->fitModel(model_type, model_params, transformations);
			}
			alignment->transformPeakMaps(peak_maps, transformations);

			// write output
			progresslogger.startProgress(0, outs.size(), "writing output files");
			for (Size i = 0; i < outs.size(); ++i)
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
			progresslogger.startProgress(0, ins.size(), "loading input files");
			for (Size i = 0; i < ins.size(); ++i)
			{
				progresslogger.setProgress(i);
		    f.load(ins[i], feat_maps[i]);
			}
			progresslogger.endProgress();

			// try to align
			try
			{
				alignment->alignFeatureMaps(feat_maps, transformations);
			}
			catch (Exception::NotImplemented&)
			{
				writeLog_("Error: The algorithm '" + alignment->getName() + "' cannot be used for feature data!");
				return INTERNAL_ERROR;
			}
			if (model_type != "none")
			{
				alignment->fitModel(model_type, model_params, transformations);
			}
			alignment->transformFeatureMaps(feat_maps, transformations);

			// write output
			progresslogger.startProgress(0, outs.size(), "writing output files");
			for (Size i = 0; i < outs.size(); ++i)
			{
				progresslogger.setProgress(i);
				
				//annotate output with data processing info
				addDataProcessing_(feat_maps[i], getProcessingInfo_(DataProcessing::ALIGNMENT));

		    f.store(outs[i], feat_maps[i]);
			}
			progresslogger.endProgress();
		}
    //-------------------------------------------------------------
    // perform consensus alignment
    //-------------------------------------------------------------
		else if (in_type == FileTypes::CONSENSUSXML)
		{
			// load input
			std::vector<ConsensusMap> cons_maps(ins.size());
			ConsensusXMLFile f;
			// f.setLogType(log_type_); // TODO
			progresslogger.startProgress(0, ins.size(), "loading input files");
			for (Size i = 0; i < ins.size(); ++i)
			{
				progresslogger.setProgress(i);
		    f.load(ins[i], cons_maps[i]);
			}
			progresslogger.endProgress();

			// try to align
			try
			{
				alignment->alignConsensusMaps(cons_maps, transformations);
			}
			catch (Exception::NotImplemented&)
			{
				writeLog_("Error: The algorithm '" + alignment->getName() + "' cannot be used for consensus feature data!");
				return INTERNAL_ERROR;
			}
			if (model_type != "none")
			{
				alignment->fitModel(model_type, model_params, transformations);
			}
			alignment->transformConsensusMaps(cons_maps, transformations);

			// write output
			progresslogger.startProgress(0, outs.size(), "writing output files");
			for (Size i = 0; i < outs.size(); ++i)
			{
				progresslogger.setProgress(i);
				
				//annotate output with data processing info
				addDataProcessing_(cons_maps[i], getProcessingInfo_(DataProcessing::ALIGNMENT));

		    f.store(outs[i], cons_maps[i]);
			}
			progresslogger.endProgress();
		}
    //-------------------------------------------------------------
    // perform peptide alignment
    //-------------------------------------------------------------
		else if (in_type == FileTypes::IDXML)
		{
			// load input
			std::vector< std::vector<ProteinIdentification> > protein_ids_vec(ins.size());
			std::vector< std::vector<PeptideIdentification> > peptide_ids_vec(ins.size());

			IdXMLFile f;
			// f.setLogType_(log_type_);

			progresslogger.startProgress(0, ins.size(), "loading input files");
			for (Size i = 0; i < ins.size(); ++i)
			{
				progresslogger.setProgress(i);
		    f.load(ins[i], protein_ids_vec[i], peptide_ids_vec[i]);
			}
			progresslogger.endProgress();

			// try to align
			try
			{
				alignment->alignPeptideIdentifications(peptide_ids_vec, transformations);
			}
			catch (Exception::NotImplemented&)
			{
				writeLog_("Error: The algorithm '" + alignment->getName() + "' cannot be used for peptide data!");
				return INTERNAL_ERROR;
			}
			if (model_type != "none")
			{
				alignment->fitModel(model_type, model_params, transformations);
			}
			alignment->transformPeptideIdentifications(peptide_ids_vec, 
																								 transformations);

			// write output
			progresslogger.startProgress(0, outs.size(), "writing output files");
			for (Size i = 0; i < outs.size(); ++i)
			{
				progresslogger.setProgress(i);
		    f.store(outs[i], protein_ids_vec[i], peptide_ids_vec[i]);
			}
			progresslogger.endProgress();
		}
		else
		{
			// TODO can this really happen? I think it is tested above. Otherwise
			// throw an appropriate exception?
			return ILLEGAL_PARAMETERS;
		}

		if (!trafos.empty())
		{
			for (Size i = 0; i < transformations.size(); ++i)
			{
				TransformationXMLFile().store(trafos[i], transformations[i]);
			}
		}

		return EXECUTION_OK;
	}
};

/// @endcond

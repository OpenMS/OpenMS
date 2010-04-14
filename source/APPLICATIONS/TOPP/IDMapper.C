// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_IDMapper IDMapper

	Assigns protein/peptide identifications to features or consensus features.

	This tool is typically used before @ref TOPP_ConsensusID.

	The mapping is based on retention times and mass-to-charge values. Roughly, a peptide identification is assigned to a (consensus) feature if its position lies within the boundaries of the feature or close enough to the feature centroid.
	Peptide identifications that don't match anywhere are still recorded in the resulting map, as "unassigned peptides". Protein identifications are annotated to the whole map, i.e. not to any particular (consensus) feature.

	On the peptide side, two sources for m/z values are possible (see parameter @p mz_reference): 1. m/z of the precursor of the MS2 spectrum that gave rise to the peptide identification; 2. theoretical masses computed from the amino acid sequences of peptide hits.
	(When using theoretical masses, make sure that peptide modifications were identified correctly. OpenMS currently "forgets" mass shifts that it can't assign to modifications - if that happens, masses computed from peptide sequences will be off.)

	In all cases, tolerance in RT and m/z dimension is applied according to the parameters @p rt_tolerance and @p mz_tolerance. Tolerance is understood as "plus or minus x", so the matching range is actually increased by twice the tolerance value.
	
	If several features or consensus features overlap the position of a peptide identification (taking the allowed tolerances into account), the identification is annotated to all of them.

	<B>Annotation of feature maps (featureXML input):</B>\n
	If @em all features have at least one convex hull, peptide positions are matched against the bounding boxes of the convex hulls (of individual mass traces, if available) by default. If not, the positions of the feature centroids are used. The respective coordinates of the centroids are also used for matching (in place of the corresponding ranges from the bounding boxes) if @p use_centroid_rt or @p use_centroid_mz are true.

	<B>Annotation of consensus maps (consensusXML input):</B>\n
	Peptide positions are always matched against centroid positions. By default, the consensus centroids are used. However, if @p use_subelements is set, the centroids of sub-features are considered instead. In this case, a peptide identification is mapped to a consensus feature if any of its sub-features matches.

	@deprecated The parameter handling of this tool has been reworked. For greater consistency with other tools, the parameters @p rt_delta and @p mz_delta have been renamed to @p rt_tolerance and @p mz_tolerance. The possible values of the @p mz_reference parameter have also been renamed. The default value of @p mz_tolerance has been increased from 1 ppm to a more realistic 20 ppm.\n
	Most importantly, the @p use_centroids parameter from previous versions has been split into two parameters, @p use_centroid_rt and @p use_centroid_mz. In OpenMS 1.6, peptide identifications would be matched only against monoisotopic mass traces of features if @p mz_reference was @p PeptideMass; otherwise, all mass traces would be used. This implicit behaviour has been abandoned, you can now explicitly control it with the @p use_centroid_mz parameter. @p use_centroid_mz does not take into account m/z deviations in the monoisotopic mass trace, but this can be compensated by increasing @p mz_tolerance. The new implementation should work correctly even if the monoisotopic mass trace itself was not detected.


	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_IDMapper.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDMapper
	: public TOPPBase
{

	public:

		TOPPIDMapper()
			: TOPPBase("IDMapper", "Assigns protein/peptide identifications to features or consensus features.")
		{
		}

	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("id", "<file>", "", "Protein/peptide identifications file");
			setValidFormats_("id",StringList::create("idXML"));
			registerInputFile_("in", "<file>", "", "Feature map/consensus map file");
			setValidFormats_("in",StringList::create("featureXML,consensusXML"));
			registerOutputFile_("out", "<file>", "", "Output file (the format depends on the input file format).");
			setValidFormats_("out",StringList::create("featureXML,consensusXML"));

			addEmptyLine_();
			IDMapper mapper;
			Param p = mapper.getParameters();
			registerDoubleOption_("rt_tolerance", "<value>", p.getValue("rt_tolerance"), "RT tolerance (in seconds) for the matching of peptide identifications and (consensus) features.\nTolerance is understood as 'plus or minus x', so the matching range increases by twice the given value.", false);
			setMinFloat_("rt_tolerance", 0.0);
			registerDoubleOption_("mz_tolerance", "<value>", p.getValue("mz_tolerance"), "m/z tolerance (in ppm or Da) for the matching of peptide identifications and (consensus) features.\nTolerance is understood as 'plus or minus x', so the matching range increases by twice the given value.", false);
			setMinFloat_("mz_tolerance", 0.0);
			registerStringOption_("mz_measure", "<choice>", p.getEntry("mz_measure").valid_strings[0], "Unit of 'mz_tolerance'.", false);
			setValidStrings_("mz_measure", p.getEntry("mz_measure").valid_strings);
			registerStringOption_("mz_reference","<choice>", p.getEntry("mz_reference").valid_strings[0], "Source of m/z values for peptide identifications. If 'precursor', the precursor-m/z from the idXML is used. If 'peptide',\nmasses are computed from the sequences of peptide hits; in this case, an identification matches if any of its hits matches.\n('peptide' should be used together with 'use_centroid_mz' to avoid false-positive matches.)", false);
			setValidStrings_("mz_reference", p.getEntry("mz_reference").valid_strings);
			addEmptyLine_();
			addText_("Additional options for featureXML input:");
			registerFlag_("use_centroid_rt", "Use the RT coordinates of the feature centroids for matching, instead of the RT ranges of the features/mass traces.");
			registerFlag_("use_centroid_mz", "Use the m/z coordinates of the feature centroids for matching, instead of the m/z ranges of the features/mass traces.\n(If you choose 'peptide' as 'mz_reference', you should usually set this flag to avoid false-positive matches.)");

			addEmptyLine_();
			addText_("Additional options for consensusXML input:");
			registerFlag_("use_subelements", "Match using RT and m/z of sub-features instead of consensus RT and m/z. A consensus feature matches if any of its sub-features matches.");
		}

		ExitCodes main_(int , const char**)
		{
			// LOG_DEBUG << "Starting..." << endl;
			String in = getStringOption_("in");
			FileTypes::Type in_type = FileHandler::getType(in);
			String out = getStringOption_("out");

			//----------------------------------------------------------------
			// load idXML
			//----------------------------------------------------------------
			// LOG_DEBUG << "Loading idXML..." << endl;
			vector<ProteinIdentification> protein_ids;
			vector<PeptideIdentification> peptide_ids;
			String document_id;
			IdXMLFile().load(getStringOption_("id"),protein_ids,peptide_ids, document_id);

			//----------------------------------------------------------------
			//create mapper
			//----------------------------------------------------------------
			// LOG_DEBUG << "Creating mapper..." << endl;
			IDMapper mapper;
			Param p = mapper.getParameters();
			p.setValue("rt_tolerance", getDoubleOption_("rt_tolerance"));
			p.setValue("mz_tolerance", getDoubleOption_("mz_tolerance"));
			p.setValue("mz_measure", getStringOption_("mz_measure"));
			p.setValue("mz_reference", getStringOption_("mz_reference"));
			mapper.setParameters(p);

			//----------------------------------------------------------------
			// consensusXML
			//----------------------------------------------------------------
			if (in_type == FileTypes::CONSENSUSXML)
			{
				// LOG_DEBUG << "Processing consensus map..." << endl;
				ConsensusXMLFile file;
				ConsensusMap map;
				file.load(in,map);
				
				bool measure_from_subelements=getFlag_("use_subelements");
				
				mapper.annotate(map,peptide_ids,protein_ids,measure_from_subelements);
				
				//annotate output with data processing info
				addDataProcessing_(map, getProcessingInfo_(DataProcessing::IDENTIFICATION_MAPPING));
				
				file.store(out,map);
			}

			//----------------------------------------------------------------
			// featureXML
			//----------------------------------------------------------------
			if (in_type == FileTypes::FEATUREXML)
			{
				// LOG_DEBUG << "Processing feature map..." << endl;
				FeatureMap<> map;
				FeatureXMLFile file;
				file.load(in, map);

				mapper.annotate(map, peptide_ids, protein_ids, 
												getFlag_("use_centroid_rt"), 
												getFlag_("use_centroid_mz"));
				
				//annotate output with data processing info
				addDataProcessing_(map, getProcessingInfo_(DataProcessing::IDENTIFICATION_MAPPING));
				
				file.store(out,map);
			}

			// LOG_DEBUG << "Done." << endl;
			return EXECUTION_OK;
		}

};


int main( int argc, const char** argv )
{
	TOPPIDMapper tool;
	return tool.main(argc, argv);
}

/// @endcond

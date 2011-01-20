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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SeedListGenerator.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_SeedListGenerator SeedListGenerator
	
	@brief Application to generate seed lists for feature detection.
	
<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ SeedListGenerator \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=3> @ref TOPP_FeatureFinder @n ("centroided" algorithm) </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinker </td>
		</tr>
	</table>
</CENTER>

	In feature detection algorithms, an early step is generally to identify points of interest in the LC-MS map (so-called seeds) that may later be extended to features. If supported by the feature detection algorithm (currently only the "centroided" algorithm), user-supplied seed lists allow greater control over this process.

	The SeedListGenerator can automatically create seed lists from a variety of sources. The lists are exported in featureXML format - suitable as input to FeatureFinder -, but can be converted to or from text formats using the @ref TOPP_TextExporter (with "-minimal" option to convert to CSV) and @ref TOPP_FileConverter (to convert from CSV) tools.

	What are possible use cases for custom seed lists?
	- In analyses that can take into account only features with peptide annotations, it may be useful to focus directly on certain locations in the LC-MS map - on all MS2 precursors (mzML input), or on precursors whose fragment spectra could be matched to a peptide sequence (idXML input).
	- When additional information becomes available during an analysis, one might want to perform a second, targeted round of feature detection on the experimental data. For example, once a feature map is annotated with peptide identifications, it is possible to go back to the LC-MS map and look for features near unassigned peptides, potentially with a lower score threshold (featureXML input).
	- Similarly, when features from different experiments are aligned and grouped, the consensus map may reveal where features were missed in the initial detection round in some experiments. The locations of these "holes" in the consensus map can be compiled into seed lists for the individual experiments (consensusXML input). (Note that the resulting seed lists use the retention time scale of the consensus map, which might be different from the original time scales of the experiments if e.g. the @ref TOPP_MapAligner tool was used to perform retention time correction as part of the alignment process. In this case, the RT transformations from the alignment must be applied to the LC-MS maps prior to the seed list-based feature detection runs.)


	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_SeedListGenerator.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{

  class TOPPSeedListGenerator 
  	: public TOPPBase
  {
	public:
		TOPPSeedListGenerator() :
			TOPPBase("SeedListGenerator", "Generates seed lists for feature detection.")
      {
      }

	protected:
	
		void registerOptionsAndFlags_()
      {
				registerInputFile_("in", "<file>", "",
													 "Input file (see below for details)");
				setValidFormats_("in", StringList::
												 create("mzML,idXML,featureXML,consensusXML"));
        registerOutputFileList_("out", "<file(s)>", StringList(), "Output file(s)");
        setValidFormats_("out", StringList::create("featureXML"));
				addEmptyLine_();		
				addText_("Options for idXML input:");
				registerFlag_("use_peptide_mass", "Use the monoisotopic mass of the best peptide hit for the m/z position (default: use precursor m/z)");
				addEmptyLine_();
				addText_("If input is consensusXML, one output file per constituent map is required (same order as in the consensusXML);\notherwise, exactly one output file.");
				addEmptyLine_();
				addText_("Seed lists can be generated from the file types below. The seeds are created at the indicated positions (RT/MZ):");
				addText_("- mzML: locations of MS2 precursors");
				addText_("- idXML: locations of peptide identifications");
				addText_("- featureXML: locations of unassigned peptide identifications");
				addText_("- consensusXML: locations of consensus features that do not contain sub-features from the respective map");
      }

		
		ExitCodes main_(int, const char**)
      {
				String in = getStringOption_("in");
				StringList out = getStringList_("out");
				SeedListGenerator seed_gen;
				// results (actually just one result, except for consensusXML input):
				Map<UInt64, SeedListGenerator::SeedList> seed_lists;

				Size num_maps = 0;
				FileTypes::Type in_type = FileHandler::getType(in);

				if (in_type == FileTypes::CONSENSUSXML)
				{
					ConsensusMap consensus;
					ConsensusXMLFile().load(in, consensus);
					num_maps = consensus.getFileDescriptions().size();
					if (out.size() != num_maps)
					{
						writeLog_("Error: expected " + String(num_maps) +
											" output filenames");
						return ILLEGAL_PARAMETERS;
					}
					seed_gen.generateSeedLists(consensus, seed_lists);
				}

				else if (out.size() > 1)
				{
					writeLog_("Error: expected only one output filename");
					return ILLEGAL_PARAMETERS;
				}
				
				else if (in_type == FileTypes::MZML)
				{
					MSExperiment<> experiment;
					MzMLFile().load(in, experiment);
					seed_gen.generateSeedList(experiment, seed_lists[0]);
				}

				else if (in_type == FileTypes::IDXML)
				{
					vector<ProteinIdentification> proteins;
					vector<PeptideIdentification> peptides;
					IdXMLFile().load(in, proteins, peptides);
					seed_gen.generateSeedList(peptides, seed_lists[0], 
																		getFlag_("use_peptide_mass"));
				}
				
				else if (in_type == FileTypes::FEATUREXML)
				{
          FeatureMap<> features;
          FeatureXMLFile().load(in, features);
					seed_gen.generateSeedList(
						features.getUnassignedPeptideIdentifications(), seed_lists[0]);
				}

				// output:
				num_maps = 0;
				for (Map<UInt64, SeedListGenerator::SeedList>::Iterator it =
							 seed_lists.begin(); it != seed_lists.end(); ++it, ++num_maps)
				{
					FeatureMap<> features;
					seed_gen.convertSeedList(it->second, features);
					//annotate output with data processing info:
					addDataProcessing_(features, getProcessingInfo_(
															 DataProcessing::DATA_PROCESSING));
					FeatureXMLFile().store(out[num_maps], features);
				}
		
				return EXECUTION_OK;
			}
	};
	
} 


int main(int argc, const char** argv)
{
  TOPPSeedListGenerator t;
  return t.main(argc, argv);
}

/// @endcond

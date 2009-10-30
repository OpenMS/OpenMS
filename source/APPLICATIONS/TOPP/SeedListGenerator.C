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
	
	Seed lists specify locations in an MS experiment where features are expected. Currently, only the "wavelet" FeatureFinder algorithm supports custom seed lists (in FeatureXML format).

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

		void writeCSV_(const String& filename,
									 const SeedListGenerator::SeedList& seeds)
			{
				ofstream outstr(filename.c_str());
				SVOutStream output(outstr, ",");
				output.modifyStrings(false);
				output << "#rt" << "mz" << endl;
				output.modifyStrings(true);
				for (SeedListGenerator::SeedList::const_iterator it = seeds.begin();
						 it != seeds.end(); ++it)
				{
					output << it->getX() << it->getY() << endl;
				}
				outstr.close();
			}
		
		
		void registerOptionsAndFlags_()
      {
				registerInputFile_("in", "<file>", "",
													 "Input file (see below for details)");
				setValidFormats_("in", StringList::
												 create("mzML,idXML,featureXML,consensusXML"));
        registerOutputFileList_("out", "<file(s)>", StringList(), "Output file(s) in featureXML or text/CSV format (detected by file extension)");
				addEmptyLine_();
				addText_("If input is consensusXML, one output file per constituent map is required (same order as in the consensusXML).\nOtherwise, one output file is expected.");
        // setValidFormats_("out", StringList::create("featureXML"));
				addEmptyLine_();
				addText_("Seed lists can be generated from the file types below; seeds are created at the indicated positions (RT/MZ):");
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

				FileTypes::Type in_type = FileHandler::getType(in);
				if (in_type == FileTypes::CONSENSUSXML)
				{
					ConsensusMap consensus;
					ConsensusXMLFile().load(in, consensus);
					Size num_maps = consensus.getFileDescriptions().size();
					if (out.size() != num_maps)
					{
						writeLog_("Error: expected " + String(num_maps) +
											" output filenames");
						return ILLEGAL_PARAMETERS;
					}

					Map<UInt64, SeedListGenerator::SeedList> seed_lists;
					seed_gen.getSeedLists(consensus, seed_lists);
					num_maps = 0;
					for (Map<UInt64, SeedListGenerator::SeedList>::Iterator it =
								 seed_lists.begin(); it != seed_lists.end(); ++it, ++num_maps)
					{
						FileTypes::Type out_type = FileHandler::getType(out[num_maps]);
						if (out_type == FileTypes::FEATUREXML)
						{
							FeatureMap<> features;
							seed_gen.convert(it->second, features);
							FeatureXMLFile().store(out[num_maps], features);
						}
						else
						{
							writeCSV_(out[num_maps], it->second);
						}			
					}
					return EXECUTION_OK;
				}

				if (out.size() > 1)
				{
					writeLog_("Error: expected only one output filename");
					return ILLEGAL_PARAMETERS;
				}
				
				SeedListGenerator::SeedList seeds;
				
				if (in_type == FileTypes::MZML)
				{
					MSExperiment<> experiment;
					MzMLFile().load(in, experiment);
					seed_gen.getSeedList(experiment, seeds);
				}

				else if (in_type == FileTypes::IDXML)
				{
					vector<ProteinIdentification> proteins;
					vector<PeptideIdentification> peptides;
					IdXMLFile().load(in, proteins, peptides);
					seed_gen.getSeedList(peptides, seeds);
				}
				
				else if (in_type == FileTypes::FEATUREXML)
				{
          FeatureMap<> features;
          FeatureXMLFile().load(in, features);
					seed_gen.getSeedList(features.getUnassignedPeptideIdentifications(),
															 seeds);
				}
				
				FileTypes::Type out_type = FileHandler::getType(out[0]);
				if (out_type == FileTypes::FEATUREXML)
				{
					FeatureMap<> features;
					seed_gen.convert(seeds, features);
					FeatureXMLFile().store(out[0], features);
				}
				else
				{
					writeCSV_(out[0], seeds);
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

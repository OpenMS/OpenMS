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

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>

#include <gsl/gsl_statistics.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page IDMassAccuracy IDMassAccuracy
	
	@brief This small utility can create decoy databases used for decoy database searches.
		
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_IDMassAccuracy.cli

	Given a number of peak maps and for each of the maps an idXML file which contains
	peptide identifications the theoretical masses of the identifications and the peaks
	of the spectra are compared. This can be done for precursor information stored in
	the spectra as well as for fragment information.

	The result is a distribution of errors of experimental vs. theoretical masses, which 
	can be used for visualization with HistView for example. Having such distributions given
	the search parameters of the sequence database search can be adjusted to speed-up 
	the identification process and to get a higher performance.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

struct MassDifference
{
	DoubleReal exp_mz;
	Int charge;
	DoubleReal theo_mz;
	DoubleReal theo_mass;
	DoubleReal intensity;
};

class TOPPIDMassAccuracy
	: public TOPPBase
{
	public:
		TOPPIDMassAccuracy()
			: TOPPBase("IDMassAccuracy","Given mass spectra and IDs a distribution of the mass error is written.", false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFileList_("in","<file list>", StringList(), "Input mzML file list, containing the spectra.");
			registerInputFileList_("id_in", "<file list>", StringList(), "Input idXML file list, containing the identifications.");

			registerOutputFile_("precursor_out","<file>","","Output file which contains the deviations from the precursors", false, false);
			registerStringList_("precursor_columns", "<columns>", StringList::create("MassDifference"), "Columns which will be written to the output file");
			setValidStrings_("precursor_columns", StringList::create("MassDifference"));
			
			registerOutputFile_("fragment_out", "<file>", "", "Output file which contains the fragment ion m/z deviations", false, false);
			registerStringList_("fragment_columns", "<columns>", StringList::create("MassDifference"), "Columns which will be written to the output file");
			setValidStrings_("fragment_columns", StringList::create("MassDifference"));
			
			registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.5, "Maximal fragment mass tolerance which is allowed for MS/MS spectra, used for the calculation of matching ions.", false, false);
			registerStringOption_("separator", "<character>", "	", "character which should be used to separate the columns in the output files");
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
		
			StringList id_in(getStringList_("id_in"));
			StringList in(getStringList_("in"));

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

		  vector<vector<PeptideIdentification> > pep_ids;
		  vector<vector<ProteinIdentification> > prot_ids;
  		pep_ids.resize(id_in.size());
  		prot_ids.resize(id_in.size());

  		IdXMLFile idxmlfile;
  		for (Size i = 0; i != id_in.size(); ++i)
  		{
    		String doc_id;
    		idxmlfile.load(id_in[i], prot_ids[i], pep_ids[i], doc_id);
  		}

  		// read mzML files
  		vector<RichPeakMap> maps;
  		maps.resize(in.size());

  		if (in.size() != id_in.size())
  		{
    		writeLog_("Number of spectrum files and identification files differs...");
    		return ILLEGAL_PARAMETERS;
  		}

  		MzMLFile mzml_file;
  		for (Size i = 0; i != in.size(); ++i)
  		{
    		mzml_file.load(in[i], maps[i]);
  		}

			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------					

			// mapping ids
      IDMapper mapper;
      for (Size i = 0; i != maps.size(); ++i)
      {
        mapper.annotate(maps[i], pep_ids[i], prot_ids[i]);
      }

			// generate precursor statistics
			vector<MassDifference> precursor_diffs;
			for (Size i = 0; i != maps.size(); ++i)
			{
				for (Size j = 0; j != maps[i].size(); ++j)
				{
					if (maps[i][j].getPeptideIdentifications().size() == 0)
					{
						continue;
					}
					for (vector<PeptideIdentification>::const_iterator it = maps[i][j].getPeptideIdentifications().begin(); it != maps[i][j].getPeptideIdentifications().end(); ++it)
					{
						if (it->getHits().size() > 0)
						{
							PeptideHit hit = *it->getHits().begin();
							if (!hit.getSequence().isValid())
							{
								continue;
							}
							MassDifference md;
							Int charge = hit.getCharge();
							if (charge == 0)
							{
								charge = 1;
							}
							md.exp_mz = (double)it->getMetaValue("MZ");
							md.theo_mz = (hit.getSequence().getMonoWeight() + (double)charge)/(double)charge;
							md.theo_mass = hit.getSequence().getMonoWeight();
							md.charge = charge;
							precursor_diffs.push_back(md);
						}
					}
				}
			}


			// generate fragment ions statistics
			vector<MassDifference> fragment_diffs;
			TheoreticalSpectrumGenerator tsg;
			SpectrumAlignment sa;
			double fragment_mass_tolerance(getDoubleOption_("fragment_mass_tolerance"));
			Param sa_param(sa.getParameters());
			sa_param.setValue("tolerance", fragment_mass_tolerance);
			sa.setParameters(sa_param);
			for (Size i = 0; i != maps.size(); ++i)
			{
				for (Size j = 0; j != maps[i].size(); ++j)
				{
					if (maps[i][j].getPeptideIdentifications().size() == 0)
					{
						continue;
					}
					for (vector<PeptideIdentification>::const_iterator it = maps[i][j].getPeptideIdentifications().begin(); it != maps[i][j].getPeptideIdentifications().end(); ++it)
					{
						if (it->getHits().size() > 0)
						{
							PeptideHit hit = *it->getHits().begin();

							if (!hit.getSequence().isValid())
							{
								continue;
							}
							RichPeakSpectrum theo_spec;
							tsg.addPeaks(theo_spec, hit.getSequence(), Residue::YIon);
							tsg.addPeaks(theo_spec, hit.getSequence(), Residue::BIon);
						
							vector<pair<Size, Size> > pairs;
							sa.getSpectrumAlignment(pairs, theo_spec, maps[i][j]);
							for (vector<pair<Size, Size> >::const_iterator pit = pairs.begin(); pit != pairs.end(); ++pit)
							{
								MassDifference md;
								md.exp_mz = maps[i][j][pit->second].getMZ();
								md.theo_mz = theo_spec[pit->first].getMZ();
								md.intensity = maps[i][j][pit->second].getIntensity();
								md.charge = hit.getCharge();
								if (md.charge != 0)
								{
									md.theo_mass = (md.theo_mz * md.charge) - md.charge;
								}
								fragment_diffs.push_back(md);
							}
						}
					}
				}
			}

			
			//-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

			String precursor_out_file(getStringOption_("precursor_out"));
			if (precursor_out_file != "")
			{
				vector<DoubleReal> errors;
				ofstream precursor_out(precursor_out_file.c_str());
				for (Size i = 0; i != precursor_diffs.size(); ++i)
				{
					precursor_out << precursor_diffs[i].exp_mz - precursor_diffs[i].theo_mz << endl;
					errors.push_back(precursor_diffs[i].exp_mz - precursor_diffs[i].theo_mz);
				}
				precursor_out.close();

				cout << "Precursor mean error: " << gsl_stats_mean(&errors.front(), 1, errors.size()) << endl;
				cout << "Precursor abs. dev.:  " << gsl_stats_absdev(&errors.front(), 1, errors.size()) << endl;
				cout << "Precursor std. dev.:  " << gsl_stats_sd(&errors.front(), 1, errors.size()) << endl;
			}

			String fragment_out_file(getStringOption_("fragment_out"));
			if (fragment_out_file != "")
			{
				vector<DoubleReal> errors;
				ofstream fragment_out(fragment_out_file.c_str());
				for (Size i = 0; i != fragment_diffs.size(); ++i)
				{
					fragment_out << fragment_diffs[i].exp_mz - fragment_diffs[i].theo_mz << endl;
					errors.push_back(fragment_diffs[i].exp_mz - fragment_diffs[i].theo_mz);
				}
				fragment_out.close();
				cout << "Fragment mean error:  " << gsl_stats_mean(&errors.front(), 1, errors.size()) << endl;
				cout << "Fragment abs. dev.:   " << gsl_stats_absdev(&errors.front(), 1, errors.size()) << endl;
				cout << "Fragment std. dev.:   " << gsl_stats_sd(&errors.front(), 1, errors.size()) << endl;
			}
	
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPIDMassAccuracy tool;
	return tool.main(argc,argv);
}
  
/// @endcond






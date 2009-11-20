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

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <algorithm>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_PrecursorMassCorrector PrecursorMassCorrector
	
	@brief Corrects the precursor entries of MS/MS spectra, by using MS1 information.

	There is currently only one option implemented. The tool expects as input data produced from
	an ABI/SCIEX QStar Pulsar I, or a similar instrument. It is a QTof instrument with 
	a moderate resolution. 

	@experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PrecursorMassCorrector.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPrecursorMassCorrector
	: public TOPPBase
{
	public:
		TOPPPrecursorMassCorrector()
			: TOPPBase("PrecursorMassCorrector","Corrects the precursor entries of MS/MS spectra, by using MS1 information.", false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","Input mzML file containing the spectra.");
			setValidFormats_("in", StringList::create("mzML"));
			registerOutputFile_("out","<file>","","Output mzML file.");
			setValidFormats_("in", StringList::create("mzML"));
			registerStringOption_("type", "<instrument_type>", "", "Type of instrument used");
			setValidStrings_("type", StringList::create("QStarPulsarI"));
			registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 1.5, "Maximal deviation in Th which is acceptable to be corrected.", false);
			setMinFloat_("precursor_mass_tolerance", 0);

			registerIntOption_("max_charge", "<charge>", 3, "Maximal charge that should be assumend for precursor peaks", false, true);
			registerDoubleOption_("intensity_threshold", "<threshold>", -1.0, "Intensity threshold value for isotope wavelet feature finder, please look at the documentation of the class for details.", false, true);
			/// @improvement add merging of spectra of the same precursor and feature

		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
			DoubleReal precursor_mass_tolerance(getDoubleOption_("precursor_mass_tolerance"));
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			PeakMap exp;
			MzMLFile().load(in, exp);

			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------					

			FeatureFinderAlgorithmIsotopeWavelet<Peak1D, Feature> iso_ff;
      Param ff_param(iso_ff.getParameters());
      ff_param.setValue("max_charge", getIntOption_("max_charge"));
      ff_param.setValue("intensity_threshold", getDoubleOption_("intensity_threshold"));
      iso_ff.setParameters(ff_param);

      FeatureFinder ff;
      ff.setLogType(ProgressLogger::NONE);

			PeakMap exp2 = exp;
			exp2.clear(false);
			for (PeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it)
			{
				if (it->size() != 0)
				{
					exp2.push_back(*it);
				}
			}

			exp = exp2;
			exp.updateRanges();

			// TODO check MS2 and MS1 counts
			ProgressLogger progresslogger;
			progresslogger.setLogType(log_type_);
			progresslogger.startProgress(0, exp.size(), "Correcting precursor masses");
			Size counter(0);
			for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
			{
				progresslogger.setProgress(++counter);
				if (it->getMSLevel() == 2)
				{
					// find first MS1 scan of the MS/MS scan
					PeakMap::Iterator ms1_it = it;
					while (ms1_it != exp.begin() && ms1_it->getMSLevel() != 1)
					{
						--ms1_it;
					}
					if (ms1_it == exp.begin() && ms1_it->getMSLevel() != 1)
					{
						writeLog_("Did not find a MS1 scan to the MS/MS scan at RT=" + String(it->getRT()));
						continue;
					}
					if (ms1_it->size() == 0)
					{
						writeDebug_("No peaks in scan at RT=" + String(ms1_it->getRT()) + String(", skipping"), 1);
						continue;
					}

					PeakMap new_exp;
					new_exp.push_back(*ms1_it);
					new_exp.updateRanges();
          FeatureMap<> features, seeds;
          ff.run("isotope_wavelet", new_exp, features, ff_param, seeds);

					if (features.size() == 0)
					{
						writeDebug_("No features found for scan RT=" + String(ms1_it->getRT()), 1);
						continue;
					}

					PeakMap::Iterator ms2_it = ms1_it;
					++ms2_it;

					while (ms2_it != exp.end() && ms2_it->getMSLevel() == 2)
					{
						if (ms2_it->getPrecursors().size() == 0)
						{
							writeDebug_("Warning: found no precursors of spectrum RT=" + String(ms2_it->getRT()) + ", skipping it.", 1);
							++ms2_it;
							continue;
						}
						else if (ms2_it->getPrecursors().size() > 1)
						{
							writeLog_("Warning: found more than one precursor of spectrum RT=" + String(ms2_it->getRT()) + ", using first one.");
						}

						Precursor prec = *ms2_it->getPrecursors().begin();
						DoubleReal prec_pos = prec.getMZ();
						
						DoubleReal min_dist(numeric_limits<DoubleReal>::max());
						Size min_feat_idx(0);
						for (Size i = 0; i != features.size(); ++i)
						{
							if (fabs(features[i].getMZ() - prec_pos) < min_dist)
							{
								min_feat_idx = i;
								min_dist = fabs(features[i].getMZ() - prec_pos);
							}
						}

						writeDebug_(" min_dist=" + String(min_dist) + " mz=" + String(features[min_feat_idx].getMZ()) + " charge=" + String(features[min_feat_idx].getCharge()), 5);
						if (min_dist < precursor_mass_tolerance)
						{
							prec.setMZ(features[min_feat_idx].getMZ());
							prec.setCharge(features[min_feat_idx].getCharge());
							vector<Precursor> precs;
							precs.push_back(prec);
							ms2_it->setPrecursors(precs);
							writeDebug_("Correcting precursor mass of spectrum RT=" + String(ms2_it->getRT()) + " from " + String(prec_pos) + " to " + String(prec.getMZ()) + " (z=" + String(prec.getCharge()) + ")", 1);
						}
						
						++ms2_it;
					}
				}
			}
			progresslogger.endProgress();

			//-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
		
			MzMLFile().store(out, exp);
	
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPPrecursorMassCorrector tool;
	return tool.main(argc,argv);
}
  
/// @endcond






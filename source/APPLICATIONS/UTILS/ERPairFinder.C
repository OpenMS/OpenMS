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
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page ERPairFinder ERPairFinder
	
	@brief Util which can be used to evaluate pair ratios on enhanced resolution (zoom) scans

	This tool allows to calculate ratios of labeled peptides based on single enhanced resolution
	scans (also called zoom scans). Zoom scans is a mode of some mass spectrometry instruments
	to allow scanning at a higher resolution, at the cost of low scan speed. It can be used to
	determine charge states of precursors on ion-trap or related instruments. 

	This tool works scan based. Each scan is examined using the IsotopeWavelet (see docs of that
	class) to fit an isotope distribution based on the averagine model. Known pairs are given
	by the pair_in input parameter, which allow to search for specific pairs. Light and heavy
	variant are search for, and the pairs are finally reported with their ratios.

	If a pair is available in several scans, the intensities are summed up and the ratio is 
	calculated from the sum of the isotope fits.

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_ERPairFinder.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

// simple helper struct which stores
// a SILAC pair, with m/z value rt
// charge 
struct SILAC_pair
{   
  DoubleReal mz_light;
	DoubleReal mz_heavy;
	Int charge;
	DoubleReal rt;
};  


// helper struct which stores the 
// SILAC_pair which it is matched to
struct MatchedFeature
{
	MatchedFeature(const Feature& feature, Size index)
		: f(feature),
			idx(index)
	{
	}

	Feature f;
	Size idx;
};

// This struct store quantitation for one scan
// for fast access to defined pair
struct SILACQuantitation
{
	SILACQuantitation(DoubleReal l_intensity, DoubleReal h_intensity, Size index)
		:	light_intensity(l_intensity),
			heavy_intensity(h_intensity),
			idx(index)
	{
	}
	DoubleReal light_intensity;
	DoubleReal heavy_intensity;
	Size idx;
};

class TOPPERPairFinder
	: public TOPPBase
{
	public:
		TOPPERPairFinder()
			: TOPPBase("ERPairFinder","Bla...", false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","Input mzML file containing the ER spectra.");
			setValidFormats_("in", StringList::create("mzML"));

			registerInputFile_("pair_in", "<file>", "", "Pair-file in the format: m/z-light m/z-heavy charge rt");

			registerOutputFile_("out","<file>","","Output consensusXML file were the decoy database will be written to.");
			setValidFormats_("out", StringList::create("consensusXML"));

			registerDoubleOption_("lower_mz", "<m/z>", 200.0, "Lower m/z value for the resampling, which is useful for the stability of the method.", false, true);
			setMinFloat_("lower_mz", 0.0);

			registerDoubleOption_("upper_mz", "<m/z>", 2000.0, "Upper m/z value for the resampling of the scans, which is useful for the stability of the methode.", false, true);
			setMinFloat_("upper_mz", 0.0);

			registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 0.3, "Precursor mass tolerance which is used for the pair finding and the matching of the given pair m/z values to the features.", false, false);
			setMinFloat_("precursor_mass_tolerance", 0.0);

			registerDoubleOption_("RT_tolerance", "<tolerance>", 200, "Maximal deviation in RT dimension in second a feature can have when comparing to the RT values given in the pair file", false, false);
			setMinFloat_("RT_tolerance", 1.0);

			registerIntOption_("max_charge", "<charge>", 3, "Maximal charge state features should be search for.", false, false);
			setMinInt_("max_charge", 1);

			registerDoubleOption_("intensity_threshold", "<threshold>", -1.0, "Intensity threshold, for the meaning see the documentation of the IsotopeWaveletFeatureFinder documentation.", false, false);
			setMinFloat_("intensity_threshold", -1.0);


		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
			String pair_in(getStringOption_("pair_in"));
			DoubleReal lower_mz(getDoubleOption_("lower_mz"));
			DoubleReal upper_mz(getDoubleOption_("upper_mz"));
			DoubleReal precursor_mass_tolerance(getDoubleOption_("precursor_mass_tolerance"));
			DoubleReal RT_tolerance(getDoubleOption_("RT_tolerance"));
			Int debug(getIntOption_("debug"));

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			PeakMap exp;
		  MzMLFile().load(in, exp);
		  exp.sortSpectra();
			exp.updateRanges();

  		// read pair file
  		ifstream is(pair_in.c_str());
  		String line;
  		vector<SILAC_pair> pairs;
  		while (getline(is, line))
  		{
    		line.trim();
    		if (line.size() == 0 || line[0] == '#')
    		{
      		continue;
   			}
    		vector<String> split;
    		line.split(' ', split);
    		if (split.size() != 4)
   	 		{
      		cerr << "missformated line ('" << line << "') should be (space separated) 'm/z-light m/z-heavy charge rt'" << endl;
    		}
    		SILAC_pair p;
    		p.mz_light = split[0].toDouble();
    		p.mz_heavy = split[1].toDouble();
    		p.charge = split[2].toInt();
    		p.rt = split[3].toDouble();
    		pairs.push_back(p);
  		}
			is.close();

			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------					


			ConsensusMap results_map;
			results_map.getFileDescriptions()[0].label = "light";
			results_map.getFileDescriptions()[0].filename = in;
  		results_map.getFileDescriptions()[1].label = "heavy";
			results_map.getFileDescriptions()[1].filename = in;

			FeatureFinderAlgorithmIsotopeWavelet<Peak1D, Feature> iso_ff;
			Param ff_param(iso_ff.getParameters());
			ff_param.setValue("max_charge", 3);
			ff_param.setValue("intensity_threshold", -1.0);
			iso_ff.setParameters(ff_param);
			
			LinearResampler res;
			res.setLogType(ProgressLogger::NONE);
			Param res_param(res.getParameters());
			res.setParameters(res_param);

			FeatureFinder ff;
			ff.setLogType(ProgressLogger::NONE);

			Size feature_counter(0);
			vector<SILACQuantitation> quantlets;
			FeatureMap<> all_features;
			for (PeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it)
			{
				if (it->size() == 0)
				{
					continue;
				}

				FeatureMap<> feature_map;
				PeakMap new_exp;
				PeakSpectrum new_spec = *it;

				// get spacing from data
				DoubleReal min_spacing(numeric_limits<DoubleReal>::max());
				DoubleReal last_mz(0);
				for (PeakSpectrum::ConstIterator pit = new_spec.begin(); pit != new_spec.end(); ++pit)
				{
					if (pit->getMZ() - last_mz < min_spacing)
					{
						min_spacing = pit->getMZ() - last_mz;
					}
					last_mz = pit->getMZ();
				}
				writeDebug_("Min-spacing=" + String(min_spacing), 1);
				res_param.setValue("spacing", min_spacing);

				// insert peaks at beginning and end to 
				// expand the resampling range to 'lower_mz' to 'upper_mz'
				Peak1D p;
				p.setMZ(lower_mz);
				p.setIntensity(0);
				new_spec.insert(new_spec.begin(), p);
				p.setMZ(upper_mz);
				new_spec.push_back(p);

				new_exp.push_back(new_spec);
				
				res.rasterExperiment(new_exp);
				
				writeDebug_("Spectrum-id: " + it->getNativeID() + " @ " + String(it->getRT()) +"s", 1);

				new_exp.updateRanges();
				new_exp.sortSpectra(true);

				ff.run("isotope_wavelet", new_exp, feature_map, ff_param);

				writeDebug_("#features=" + String(feature_map.size()), 1);

				// search if feature maps to m/z value of pair
				vector<MatchedFeature> light, heavy;
				for (FeatureMap<>::const_iterator fit = feature_map.begin(); fit != feature_map.end(); ++fit)
				{
					all_features.push_back(*fit);
					Size idx = 0;
					for (vector<SILAC_pair>::const_iterator pit = pairs.begin(); pit != pairs.end(); ++pit, ++idx)
					{
						if (fabs(fit->getMZ() - pit->mz_light) < precursor_mass_tolerance && fabs(fit->getRT() - pit->rt) < RT_tolerance)
						{
							light.push_back(MatchedFeature(*fit, idx));
						}
						if (fabs(fit->getMZ() - pit->mz_heavy) < precursor_mass_tolerance && fabs(fit->getRT() - pit->rt) < RT_tolerance)
						{
							heavy.push_back(MatchedFeature(*fit, idx));
						}
					}
				
				}

				if (heavy.size() != 0 && light.size() != 0)
				{
					writeDebug_("Finding best feature pair out of " + String(light.size()) + " light and " + String(heavy.size()) + " heavy matching features.", 1);
					// now find "good" matches, means the pair with the smallest m/z deviation
					Feature best_light, best_heavy;
					DoubleReal best_deviation(numeric_limits<DoubleReal>::max());
					Size best_idx(pairs.size());
					for (vector<MatchedFeature>::const_iterator fit1 = light.begin(); fit1 != light.end(); ++fit1)
					{
						for (vector<MatchedFeature>::const_iterator fit2 = heavy.begin(); fit2 != heavy.end(); ++fit2)
						{
							if (fit1->idx != fit2->idx)
							{
								continue;
							}
							DoubleReal deviation(0);
							deviation = fabs(fit1->f.getMZ() - pairs[fit1->idx].mz_light) + 
													fabs(fit2->f.getMZ() - pairs[fit2->idx].mz_heavy);
							if (deviation < best_deviation)
							{
								best_light = fit1->f;
								best_heavy = fit2->f;
								best_idx = fit1->idx;
							}
						}
					}

					if (best_idx == pairs.size())
					{
						continue;
					}

					writeDebug_("Ratio: " + String(best_heavy.getIntensity() / best_light.getIntensity()), 1);
					ConsensusFeature SILAC_feature;
					SILAC_feature.setMZ((best_light.getMZ() + best_heavy.getMZ()) / 2.0);
					SILAC_feature.setRT((best_light.getRT() + best_heavy.getRT()) / 2.0);
	    		SILAC_feature.insert(0, feature_counter, best_light);
	    		SILAC_feature.insert(1, feature_counter++, best_heavy);
  	  		results_map.push_back(SILAC_feature);

					quantlets.push_back(SILACQuantitation(best_light.getIntensity(), best_heavy.getIntensity(), best_idx));
				}
			}

			// now calculate the final quantitation values from the quantlets
			Map<Size, vector<SILACQuantitation> > idx_to_quantlet;
			for (vector<SILACQuantitation>::const_iterator it = quantlets.begin(); it != quantlets.end(); ++it)
			{
				idx_to_quantlet[it->idx].push_back(*it);
			}

			for (Map<Size, vector<SILACQuantitation> >::ConstIterator it1 = idx_to_quantlet.begin(); it1 != idx_to_quantlet.end(); ++it1)
			{
				SILAC_pair silac_pair = pairs[it1->first];
				writeDebug_("Quantitation of pair " + String(silac_pair.mz_light) + " <-> " + String(silac_pair.mz_heavy) + " @RT=" + String(silac_pair.rt) + "s (#scans for quantation=" + String(it1->second.size()) + ")", 1);

				// simply add up all intensities and calculate the final ratio
				DoubleReal light_sum(0), heavy_sum(0);
				for (vector<SILACQuantitation>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
				{
					light_sum += it2->light_intensity;
					heavy_sum += it2->heavy_intensity;
				}

				cout << "Ratio: " << silac_pair.mz_light << " <-> " << silac_pair.mz_heavy << " @ " << silac_pair.rt << " s, ratio(h/l) " << heavy_sum / light_sum << endl;
			}


			//-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

			if (debug > 1)
			{
				writeDebug_("Writing featureXML file with all the features", 2);
				FeatureXMLFile().store(out.prefix('.') + ".featureXML", all_features);
			}
			writeDebug_("Writing output", 1);
			ConsensusXMLFile().store(out, results_map);
	
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPERPairFinder tool;
	return tool.main(argc,argv);
}
  
/// @endcond






// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <gsl/gsl_statistics.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_ERPairFinder ERPairFinder
	
	@brief Util which can be used to evaluate pair ratios on enhanced resolution (zoom) scans.

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
	
	@experimental This software is experimental and might contain bugs!

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
			: TOPPBase("ERPairFinder","Util which can be used to evaluate pair ratios on enhanced resolution (zoom) scans.", false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","Input mzML file containing the ER spectra.");
			setValidFormats_("in", StringList::create("mzML"));

			registerInputFile_("pair_in", "<file>", "", "Pair-file in the format: m/z-light m/z-heavy charge rt");

			registerOutputFile_("out","<file>","","Output consensusXML file were the pairs of the feature are written into.");
			setValidFormats_("out", StringList::create("consensusXML"));

			registerOutputFile_("feature_out", "<file>", "", "Output featureXML file, only written if given, skipped otherwise.", false, false);
			setValidFormats_("feature_out", StringList::create("featureXML"));

			registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 0.3, "Precursor mass tolerance which is used for the pair finding and the matching of the given pair m/z values to the features.", false, false);
			setMinFloat_("precursor_mass_tolerance", 0.0);

			registerDoubleOption_("RT_tolerance", "<tolerance>", 200, "Maximal deviation in RT dimension in seconds a feature can have when comparing to the RT values given in the pair file", false, true);
			setMinFloat_("RT_tolerance", 1.0);

			registerIntOption_("max_charge", "<charge>", 3, "Maximal charge state features should be search for.", false, true);
			setMinInt_("max_charge", 1);

			registerDoubleOption_("intensity_threshold", "<threshold>", -1.0, "Intensity threshold, for the meaning see the documentation of the IsotopeWaveletFeatureFinder documentation.", false, true);
			setMinFloat_("intensity_threshold", -1.0);

			registerIntOption_("max_isotope", "<num>", 3, "Max isotope of the isotope distribution to be considered", false, true);
			setMinInt_("max_isotope", 2);

			registerDoubleOption_("expansion_range", "<range>", 5.0, "The range that is used to extend the isotope distribution with null intensity peaks in Th.", false, true);
			setMinFloat_("expansion_range", 0.0);
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
			String pair_in(getStringOption_("pair_in"));
			String feature_out(getStringOption_("feature_out"));
			DoubleReal precursor_mass_tolerance(getDoubleOption_("precursor_mass_tolerance"));
			DoubleReal RT_tolerance(getDoubleOption_("RT_tolerance"));
			DoubleReal expansion_range(getDoubleOption_("expansion_range"));
			Size max_isotope(getIntOption_("max_isotope"));
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
    		if (line.empty() || line[0] == '#')
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
			
			FeatureFinder ff;
			ff.setLogType(ProgressLogger::NONE);

			vector<SILACQuantitation> quantlets;
			FeatureMap<> all_features;
			for (PeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it)
			{
				if (it->size() == 0 || it->getMSLevel() != 1 || !it->getInstrumentSettings().getZoomScan())
				{
					continue;
				}

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

				// split the spectrum into two subspectra, by using different hypothesis of
				// the SILAC pairs
				Size idx = 0;
				for (vector<SILAC_pair>::const_iterator pit = pairs.begin(); pit != pairs.end(); ++pit, ++idx)
				{
					// in RT window?
					if (fabs(it->getRT() - pit->rt) >= RT_tolerance)
					{
						continue;
					}

					// now excise the two ranges for the pair, complete isotope distributions of both, light and heavy
					PeakSpectrum light_spec, heavy_spec;
					light_spec.setRT(it->getRT());
					heavy_spec.setRT(it->getRT());
					for (PeakSpectrum::ConstIterator sit = it->begin(); sit != it->end(); ++sit)
					{
						DoubleReal mz(sit->getMZ());
						if (mz - (pit->mz_light - precursor_mass_tolerance) > 0 && 
								(pit->mz_light + (DoubleReal)max_isotope * Constants::NEUTRON_MASS_U / (DoubleReal)pit->charge + precursor_mass_tolerance) - mz  > 0)
						{
							light_spec.push_back(*sit);
						}
						            
						if (mz - (pit->mz_heavy - precursor_mass_tolerance) > 0 &&      
                (pit->mz_heavy + (DoubleReal)max_isotope * Constants::NEUTRON_MASS_U / (DoubleReal)pit->charge + precursor_mass_tolerance) - mz  > 0)
            {
              heavy_spec.push_back(*sit);
            }
					}

					// expand light spectrum
        	Peak1D p;
        	p.setIntensity(0);

					if (light_spec.size() > 0)
					{
        		DoubleReal lower_border = light_spec.begin()->getMZ() - expansion_range;
        		for (DoubleReal pos = light_spec.begin()->getMZ(); pos > lower_border; pos -= min_spacing)
        		{
          		p.setMZ(pos);
          		light_spec.insert(light_spec.begin(), p);
        		}
					
        		DoubleReal upper_border = light_spec.begin()->getMZ() - expansion_range;
        		for (DoubleReal pos = light_spec.rbegin()->getMZ(); pos < upper_border; pos += min_spacing)
        		{
          		p.setMZ(pos);
          		light_spec.push_back(p);
        		}
					}

					if (heavy_spec.size() > 0)
					{
						// expand heavy spectrum
						DoubleReal lower_border = heavy_spec.begin()->getMZ() - expansion_range;
        		for (DoubleReal pos = heavy_spec.begin()->getMZ(); pos > lower_border; pos -= min_spacing)
        		{
          		p.setMZ(pos);
          		heavy_spec.insert(heavy_spec.begin(), p);
        		}

        		DoubleReal upper_border = heavy_spec.begin()->getMZ() - expansion_range;
        		for (DoubleReal pos = heavy_spec.rbegin()->getMZ(); pos < upper_border; pos += min_spacing)
        		{
          		p.setMZ(pos);
          		heavy_spec.push_back(p);
        		}
					}

					// create experiments for feature finding
					PeakMap new_exp_light, new_exp_heavy;
					new_exp_light.push_back(light_spec);
					new_exp_heavy.push_back(heavy_spec);

					if (debug > 9)
					{
						MzMLFile().store(String(it->getRT()) + "_debugging_light.mzML", new_exp_light);
						MzMLFile().store(String(it->getRT()) + "_debugging_heavy.mzML", new_exp_heavy);
					}

					writeDebug_("Spectrum-id: " + it->getNativeID() + " @ " + String(it->getRT()) +"s", 1);

					new_exp_light.updateRanges();
					new_exp_heavy.updateRanges();

					FeatureMap<> feature_map_light, feature_map_heavy, seeds;
					if (light_spec.size() > 0)
					{
						ff.run("isotope_wavelet", new_exp_light, feature_map_light, ff_param, seeds);
					}
					writeDebug_("#light_features=" + String(feature_map_light.size()), 1);
					if (heavy_spec.size() > 0)
					{
						ff.run("isotope_wavelet", new_exp_heavy, feature_map_heavy, ff_param, seeds);
					}
					writeDebug_("#heavy_features=" + String(feature_map_heavy.size()), 1);

					// search if feature maps to m/z value of pair
					vector<MatchedFeature> light, heavy;
					for (FeatureMap<>::const_iterator fit = feature_map_light.begin(); fit != feature_map_light.end(); ++fit)
					{
						all_features.push_back(*fit);
						light.push_back(MatchedFeature(*fit, idx));
					}
					for (FeatureMap<>::const_iterator fit = feature_map_heavy.begin(); fit != feature_map_heavy.end(); ++fit)
					{
						all_features.push_back(*fit);
						heavy.push_back(MatchedFeature(*fit, idx));
					}

          if ( !heavy.empty() && !light.empty() )
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
								if (fit1->idx != fit2->idx || fit1->f.getCharge() != fit2->f.getCharge() ||
										fabs(fit1->f.getMZ() - pairs[fit1->idx].mz_light) > precursor_mass_tolerance || 
										fabs(fit2->f.getMZ() - pairs[fit2->idx].mz_heavy) > precursor_mass_tolerance)
								{
									continue;
								}
								DoubleReal deviation(0);
								deviation = fabs((fit1->f.getMZ() - pairs[fit1->idx].mz_light) - (fit2->f.getMZ() - pairs[fit2->idx].mz_heavy));
								if (deviation < best_deviation && deviation < precursor_mass_tolerance)
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
	    			SILAC_feature.insert(0, best_light);
	    			SILAC_feature.insert(1, best_heavy);
  	  			results_map.push_back(SILAC_feature);
						quantlets.push_back(SILACQuantitation(best_light.getIntensity(), best_heavy.getIntensity(), best_idx));
					}
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

				// simply add up all intensities and calculate the final ratio
				DoubleReal light_sum(0), heavy_sum(0);
				vector<DoubleReal> light_ints, heavy_ints, ratios;
				for (vector<SILACQuantitation>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
				{
					light_sum += it2->light_intensity;
					light_ints.push_back(it2->light_intensity);
					heavy_sum += it2->heavy_intensity;
					heavy_ints.push_back(it2->heavy_intensity);
					ratios.push_back(it2->heavy_intensity / it2->light_intensity * (it2->heavy_intensity + it2->light_intensity));
				}

				DoubleReal absdev_ratios = gsl_stats_absdev(&ratios.front(), 1, ratios.size()) / (heavy_sum + light_sum);
				cout << "Ratio: " << silac_pair.mz_light << " <-> " << silac_pair.mz_heavy << " @ " << silac_pair.rt << " s, ratio(h/l) " << heavy_sum / light_sum << " +/- " << absdev_ratios << " (#scans for quantation: " << String(it1->second.size()) << " )" << endl;
			}


			//-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

			if (feature_out != "")
			{
				FeatureXMLFile().store(feature_out, all_features);
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






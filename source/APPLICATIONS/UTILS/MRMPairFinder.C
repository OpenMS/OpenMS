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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
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
	@page UTILS_MRMPairFinder MRMPairFinder
	
	@brief Util which can be used to evaluate pairs of MRM experiments
	
	@experimental This software is experimental and might contain bugs!
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_MRMPairFinder.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

// simple helper struct which stores
// a SILAC pair, with m/z value rt
struct SILAC_pair
{
	DoubleReal mz_light;
	DoubleReal mz_heavy;
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

class TOPPMRMPairFinder
	: public TOPPBase
{
	public:
		TOPPMRMPairFinder()
			: TOPPBase("MRMPairFinder","Util which can be used to evaluate labeled pair ratios on MRM features.", false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","Input featureXML file containing the features of the MRM experiment spectra.");
			setValidFormats_("in", StringList::create("featureXML"));

			registerInputFile_("pair_in", "<file>", "", "Pair-file in the format: prec-m/z-light prec-m/z-heavy frag-m/z-light frag-m/z-heavy rt");

			registerOutputFile_("out","<file>","","Output consensusXML file were the pairs of the features will be written to.");
			setValidFormats_("out", StringList::create("consensusXML"));

			registerOutputFile_("feature_out", "<file>", "", "Output featureXML file, only written if given, skipped otherwise.", false);
			setValidFormats_("feature_out", StringList::create("featureXML"));

			registerDoubleOption_("mass_tolerance", "<tolerance>", 0.01, "Precursor mass tolerance which is used for the pair finding and the matching of the given pair m/z values to the features.", false, true);
			setMinFloat_("mass_tolerance", 0.0);

			registerDoubleOption_("RT_tolerance", "<tolerance>", 200, "Maximal deviation in RT dimension in seconds a feature can have when comparing to the RT values given in the pair file.", false, true);
			setMinFloat_("RT_tolerance", 0.0);
			registerDoubleOption_("RT_pair_tolerance", "<tolerance>", 5, "Maximal deviation in RT dimension in seconds the two partners of a pair is allowed to have.", false, true);
			setMinFloat_("RT_pair_tolerance", 0.0);
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
			String feature_out(getStringOption_("feature_out"));
			String pair_in(getStringOption_("pair_in"));
			DoubleReal mass_tolerance(getDoubleOption_("mass_tolerance"));
			DoubleReal RT_tolerance(getDoubleOption_("RT_tolerance"));
			DoubleReal RT_pair_tolerance(getDoubleOption_("RT_pair_tolerance"));

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			FeatureMap<> all_mrm_features;
			FeatureXMLFile().load(in, all_mrm_features);

  		// read pair file
  		ifstream is(pair_in.c_str());
  		String line;
  		Map<DoubleReal, Map<DoubleReal, vector<SILAC_pair> > > pairs;
  		while (getline(is, line))
  		{
    		line.trim();
    		if (line.empty() || line[0] == '#')
    		{
      		continue;
   			}
    		vector<String> split;
    		line.split(' ', split);
				if (split.empty())
				{
					line.split('\t', split);
				}
    		if (split.size() != 5)
   	 		{
      		cerr << "missformated line ('" << line << "') should be (space separated) 'prec-m/z-light prec-m/z-heavy frag-m/z-light frag-m/z-heavy rt'" << endl;
					continue;
    		}
    		SILAC_pair p;
				DoubleReal prec_mz_light = split[0].toDouble();
				DoubleReal prec_mz_heavy = split[1].toDouble();
    		p.mz_light = split[2].toDouble();
    		p.mz_heavy = split[3].toDouble();
    		p.rt = split[4].toDouble();
    		pairs[prec_mz_light][prec_mz_heavy].push_back(p);
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

			// collect the different MRM XIC pairs for each SILAC pair as quantlets
			// then calculate the ratio over the quanlets and calculate some statistics
			FeatureMap<> all_features;
			for (Map<DoubleReal, Map<DoubleReal, vector<SILAC_pair> > >::ConstIterator it1 = pairs.begin(); it1 != pairs.end(); ++it1)
			{
				for (Map<DoubleReal, vector<SILAC_pair> >::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
				{
					vector<SILACQuantitation> quantlets;
					writeDebug_("Analyzing SILAC pair: " + String(it1->first) + " <-> " + String(it2->first), 3);
					Size idx = 0;
					for (vector<SILAC_pair>::const_iterator pit = it2->second.begin(); pit != it2->second.end(); ++pit, ++idx)
					{
						FeatureMap<> feature_map_light, feature_map_heavy;
						for (FeatureMap<>::const_iterator it = all_mrm_features.begin(); it != all_mrm_features.end(); ++it)
						{
							if (fabs((DoubleReal)it->getMetaValue("MZ") - it1->first) < mass_tolerance &&
									fabs(it->getMZ() - pit->mz_light) < mass_tolerance &&
									fabs(it->getRT() - pit->rt) < RT_tolerance)
							{
								feature_map_light.push_back(*it);
							}

							if (fabs((DoubleReal)it->getMetaValue("MZ") - it2->first) < mass_tolerance &&
              		fabs(it->getMZ() - pit->mz_heavy) < mass_tolerance &&
              		fabs(it->getRT() - pit->rt) < RT_tolerance)
          		{
            		feature_map_heavy.push_back(*it);
          		}
						}

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
							Size best_idx(it2->second.size());
							for (vector<MatchedFeature>::const_iterator fit1 = light.begin(); fit1 != light.end(); ++fit1)
							{
								for (vector<MatchedFeature>::const_iterator fit2 = heavy.begin(); fit2 != heavy.end(); ++fit2)
								{	
									if (fit1->idx != fit2->idx || fabs(fit1->f.getRT() - fit2->f.getRT()) > RT_pair_tolerance)
									{
										continue;
									}
									DoubleReal deviation(0);
									deviation = fabs(fit1->f.getMZ() - it2->second[fit1->idx].mz_light) + 
															fabs(fit2->f.getMZ() - it2->second[fit2->idx].mz_heavy);
									if (deviation < best_deviation && deviation < mass_tolerance)
									{
										best_light = fit1->f;
										best_heavy = fit2->f;
										best_idx = fit1->idx;
									}
								}
							}

							if (best_idx == it2->second.size())
							{
								continue;
							}

							ConsensusFeature SILAC_feature;
							SILAC_feature.setMZ((best_light.getMZ() + best_heavy.getMZ()) / 2.0);
							SILAC_feature.setRT((best_light.getRT() + best_heavy.getRT()) / 2.0);
	    				SILAC_feature.insert(0, best_light);
	    				SILAC_feature.insert(1, best_heavy);
  	  				results_map.push_back(SILAC_feature);

							quantlets.push_back(SILACQuantitation(best_light.getIntensity(), best_heavy.getIntensity(), best_idx));
							writeDebug_("Ratio of XIC: " + String(best_heavy.getIntensity() / best_light.getIntensity()) + " " + String(best_light.getMZ()) + " <-> " + String(best_heavy.getMZ()) + " @" + String(SILAC_feature.getRT()) + " RT-heavy=" + String(best_heavy.getRT()) + ", RT-light=" + String(best_light.getRT()) + ", RT-diff=" + String(best_heavy.getRT() - best_light.getRT()) + 
							 " avg. int " + String((best_heavy.getIntensity() + best_light.getIntensity()) / 2.0), 1);

						}
					}

					writeDebug_("Quantitation of pair " + String(it1->first) + " <-> " + String(it2->first) + " (#XIC pairs for quantation=" + String(quantlets.size()) + ")", 1);

					if (quantlets.empty())
					{
						continue;
					}

					// simply add up all intensities and calculate the final ratio
					DoubleReal light_sum(0), heavy_sum(0);
					vector<DoubleReal> light_ints, heavy_ints, ratios;
					for (vector<SILACQuantitation>::const_iterator qit1 = quantlets.begin(); qit1 != quantlets.end(); ++qit1)
					{
						light_sum += qit1->light_intensity;
						light_ints.push_back(qit1->light_intensity);
						heavy_sum += qit1->heavy_intensity;
						heavy_ints.push_back(qit1->heavy_intensity);
						ratios.push_back(qit1->heavy_intensity / qit1->light_intensity * (qit1->heavy_intensity + qit1->light_intensity));
					}

					DoubleReal absdev_ratios = gsl_stats_absdev(&ratios.front(), 1, ratios.size()) / (light_sum + heavy_sum);
					cout << "Ratio: " << it1->first << " <-> " << it2->first << " @ " << it2->second.begin()->rt << " s, ratio(h/l) " << heavy_sum / light_sum << " +/- " << absdev_ratios <<  " " << "(#XIC-pairs for quantation: " + String(ratios.size()) + " )" << endl;
				}
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
	TOPPMRMPairFinder tool;
	return tool.main(argc,argv);
}
  
/// @endcond






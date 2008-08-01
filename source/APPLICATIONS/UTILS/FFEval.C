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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FFEval FFEval
	
	@brief Evaluation tool for isotope-labeled quantitation experiments.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFFEVal
	: public TOPPBase
{
 public:
	TOPPFFEVal()
		: TOPPBase("FFEVal","Validates XML files against an XML schema.")
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		addText_("Input options");
		registerInputFile_("features","<file>","","Feature result file");
		setValidFormats_("features", StringList::create("featureXML"));
		registerInputFile_("manual","<file>","","Manual result file, which contains 6 semicolon-separated columns:\n"
																						"RT, m/z, int for feature 1 and RT, m/z, int for feature 2");
		registerInputFile_("abort_reasons","<file>","","Feature file containing the abort reasons",false);
		setValidFormats_("abort_reasons", StringList::create("featureXML"));
		registerDoubleOption_("rt_tol","<tol>",100.0,"Maximum allowed retention time deviation",false);
		registerDoubleOption_("mz_tol","<tol>",0.1,"Maximum allowed m/z deviation",false);

		addText_("Output options");
		registerOutputFile_("manual_out","<file>","","'manual' feature converted to featureXML",false);
		registerOutputFile_("unmatched_out","<file>","","unmatched 'manual' feature converted to featureXML",false);
		registerOutputFile_("ratio_out","<file>","","pairs with big ratio differences",false);
		registerOutputFile_("intensity_plot","<file>","","Gnuplot file of matched intensities",false);
		registerOutputFile_("ratio_plot","<file>","","Gnuplot file of matched intensity ratios",false);
	}

	String correlation(vector<DoubleReal> a, vector<DoubleReal> b)
	{
		DoubleReal corr = Math::pearsonCorrelationCoefficient(a.begin(),a.end(),b.begin(),b.end());
		return String::number(corr,3);
	}

	String fiveNumbers(vector<DoubleReal> a, UInt decimal_places)
	{
		sort(a.begin(),a.end());
		return String::number(a[0],decimal_places) + " " + String::number(a[a.size()/4],decimal_places) + " " + String::number(a[a.size()/2],decimal_places) + " " + String::number(a[(3*a.size())/4],decimal_places) + " " + String::number(a.back(),decimal_places);
	}
	
	String fiveNumberErrors(vector<DoubleReal> a, vector<DoubleReal> b, UInt decimal_places)
	{
		vector<DoubleReal> errors;
		for (UInt i=0; i< a.size(); ++i) errors.push_back(fabs(b[i] - a[i]));
		return fiveNumbers(errors, decimal_places);
	}

	String fiveNumberQuotients(vector<DoubleReal> a, vector<DoubleReal> b, UInt decimal_places)
	{
		vector<DoubleReal> errors;
		for (UInt i=0; i< a.size(); ++i) errors.push_back(a[i] / b[i]);
		return fiveNumbers(errors, decimal_places);
	}

	String meanError(vector<DoubleReal> a, vector<DoubleReal> b)
	{
		DoubleReal sum;
		for (UInt i=0; i< a.size(); ++i) sum += fabs(b[i] - a[i]);
		return  String::number(sum/a.size(),3);
	}

	ExitCodes main_(int , const char**)
	{
		//load automatic features
		FeatureMap<> features_automatic;
		FeatureXMLFile().load(getStringOption_("features"),features_automatic);
		features_automatic.sortByPosition();
		
		//create manual features
		TextFile text;		
		text.load(getStringOption_("manual"));
		FeatureMap<> features_manual;
		Feature tmp;		
		StringList parts;
		for(UInt i=2; i< text.size(); ++i)
		{
			text[i].split(';',parts);

			tmp.setRT(parts[0].toDouble());
			tmp.setMZ(parts[1].toDouble());		
			tmp.setIntensity(parts[2].toDouble());
			features_manual.push_back(tmp);

			tmp.setRT(parts[3].toDouble());
			tmp.setMZ(parts[4].toDouble());		
			tmp.setIntensity(parts[5].toDouble());
			features_manual.push_back(tmp);
		}
		features_manual.sortByPosition();
		if (getStringOption_("manual_out")!="")
		{
			FeatureXMLFile().store(getStringOption_("manual_out"),features_manual);
		}
		
		//load abort reasons
		FeatureMap<> abort_reasons;
		if(getStringOption_("abort_reasons")!="")
		{
			FeatureXMLFile().load(getStringOption_("abort_reasons"),abort_reasons);
		}
		
		//compare seek manual feature in automatic feature map
		DoubleReal mz_tol = getDoubleOption_("mz_tol");
		DoubleReal rt_tol = getDoubleOption_("rt_tol");
		UInt matched_single = 0;
		UInt matched_multi = 0;
		UInt matched_pairs = 0;
		vector<DoubleReal> a_int, m_int, a_ratio, m_ratio, rt_diffs, mz_diffs;
		vector<DoubleReal> matched_int, unmatched_int;
		Feature last_best_match;
		FeatureMap<> unmatched_features;
		FeatureMap<> ratio_pairs;
		Map<String,UInt> abort_strings;
		for (UInt m=0; m<features_manual.size(); ++m)
		{
			const Feature& f_m =  features_manual[m];
			//cout << "manual: " << f_m.getRT() << "/" << f_m.getMZ() << endl;
			UInt match_count = 0;
			DoubleReal best_score = 0;
			Feature best_match;
			for (UInt a=0; a<features_automatic.size(); ++a)
			{
				const Feature& f_a =  features_automatic[a];
				if ( f_a.getRT()< (f_m.getRT() + rt_tol) 
					&& f_a.getRT()> (f_m.getRT() - rt_tol)
					&& f_a.getMZ()< (f_m.getMZ() + mz_tol)
					&& f_a.getMZ()> (f_m.getMZ() - mz_tol))
				{
					++match_count;
					DoubleReal score = (1.0-fabs(f_a.getMZ()-f_m.getMZ())/mz_tol) * (1.0-fabs(f_a.getRT()-f_m.getRT())/rt_tol);
					if (score > best_score)
					{
						best_score = score;
						best_match = f_a;
					}
					//cout << " - match: " << f_a.getRT() << "/" << f_a.getMZ() << " score: " <<  score << endl;	
				}
			}
			if (match_count==1) ++matched_single;
			if (match_count>1) ++matched_multi;
			
			//calculate statistics with the best match only
			if (match_count>=1)
			{
				a_int.push_back(best_match.getIntensity());
				m_int.push_back(f_m.getIntensity());
				matched_int.push_back(f_m.getIntensity());
				//found pair
				if (Math::isOdd(m) && last_best_match!=Feature())
				{
					++matched_pairs;
					DoubleReal a_r = best_match.getIntensity()/last_best_match.getIntensity();
					a_ratio.push_back(a_r);
					DoubleReal m_r = f_m.getIntensity()/features_manual[m-1].getIntensity();
					m_ratio.push_back(m_r);
					rt_diffs.push_back(best_match.getRT()-last_best_match.getRT());
					mz_diffs.push_back(best_match.getMZ()-last_best_match.getMZ());
					
					//store pairs with too big deviation from the manual ratio
					if (a_r/m_r < 0.5 || a_r/m_r > 2.0)
					{
						ratio_pairs.push_back(f_m);
						ratio_pairs.back().setMetaValue("automatic_int", best_match.getIntensity());
						ratio_pairs.back().setMetaValue("ratio_quotient", a_r/m_r);
						ratio_pairs.back().setMetaValue("automatic_label", String(best_match.getMetaValue("label")));
						ratio_pairs.push_back(features_manual[m-1]);
						ratio_pairs.back().setMetaValue("automatic_int", last_best_match.getIntensity());
						ratio_pairs.back().setMetaValue("ratio_quotient", a_r/m_r);
						ratio_pairs.back().setMetaValue("automatic_label", String(last_best_match.getMetaValue("label")));
					}
				}
			}
			else
			{
				unmatched_int.push_back(f_m.getIntensity());
				
				//look up the abort reason of the nearest seed
				DoubleReal best_score_ab = 0;
				String reason = "";
				for (UInt b=0; b<abort_reasons.size(); ++b)
				{
					const Feature& f_ab =  abort_reasons[b];
					if ( f_ab.getRT()< (f_m.getRT() + rt_tol) 
						&& f_ab.getRT()> (f_m.getRT() - rt_tol)
						&& f_ab.getMZ()< (f_m.getMZ() + mz_tol)
						&& f_ab.getMZ()> (f_m.getMZ() - mz_tol))
					{
						DoubleReal score = (1.0-fabs(f_ab.getMZ()-f_m.getMZ())/mz_tol) * (1.0-fabs(f_ab.getRT()-f_m.getRT())/rt_tol);
						if (score > best_score_ab)
						{
							best_score_ab = score;
							reason = f_ab.getMetaValue("abort_reason");
						}
					}
				}
				if (reason=="")
				{
					reason = "No seed found";
					unmatched_features.push_back(features_manual[m]);
				}
				
				if (abort_strings.has(reason))
				{
					abort_strings[reason]++;	
				}
				else
				{
					abort_strings[reason] = 1;
				}
			}
			
			//store last best match (for pairs)
			last_best_match = best_match;
		}
		
		cout << endl;
		cout << "feature detection statistics:" << endl;
		cout << "=============================" << endl;
		cout << "  manual features: " << features_manual.size() << endl;
		cout << "  matches: " << matched_single + matched_multi << " (" << String::number(100.0*(matched_single + matched_multi)/features_manual.size(),2) << "%)" << endl;
		cout << "    one match: " << matched_single << " (" << String::number(100.0*matched_single/features_manual.size(),2) << "%)" << endl;
		cout << "    multiple matches: " << matched_multi << " (" << String::number(100.0*matched_multi/features_manual.size(),2) << "%)" << endl;
		cout << "  intensity correlation: " << correlation(m_int,a_int) << endl;
		cout << endl;
		cout << "  intensity of matched: " << fiveNumbers(matched_int,1) << endl;
		cout << "  intensity of unmatched: " << fiveNumbers(unmatched_int,1) << endl;
		cout << endl;		
		cout << "  reasons for unmatched features:" << endl;
		for (Map<String,UInt>::iterator it=abort_strings.begin(); it!=abort_strings.end(); ++it)
		{
			cout << "    " << String(it->second).fillLeft(' ',3) << ": " << it->first << endl;
		}
		cout << endl << endl;
		cout << "pair detection statistics:" << endl;
		cout << "==========================" << endl;
		cout << "  manual pairs: " <<0.5*features_manual.size() << endl;
		cout << "  found: " << matched_pairs << " (" << String::number(100.0*(matched_pairs)/(0.5*features_manual.size()),2) << "%)" << endl;
		cout << "  relative pair ratios: " << fiveNumberQuotients(m_ratio,a_ratio,3) << endl;
		cout << "  pair distance RT : " << fiveNumbers(rt_diffs, 2) << endl;
		cout << "  pair distance m/z: " << fiveNumbers(mz_diffs, 2) << endl;
		
		//write intensity pair file
		if (getStringOption_("intensity_plot")!="")
		{
			TextFile int_out;
			int_out.push_back("#manual automatic");		
			for (UInt i=0; i< a_int.size(); ++i)
			{
				int_out.push_back(String(m_int[i]) + " " + a_int[i]);		
			}
			int_out.store(getStringOption_("intensity_plot"));
		}

		//write intensity pair file
		if (getStringOption_("ratio_plot")!="")
		{
			TextFile ratio_plot;
			ratio_plot.push_back("#manual automatic");		
			for (UInt i=0; i< a_ratio.size(); ++i)
			{
				ratio_plot.push_back(String(m_ratio[i]) + " " + a_ratio[i]);		
			}
			ratio_plot.store(getStringOption_("ratio_plot"));
		}

		//write unmatched manual features
		if (getStringOption_("unmatched_out")!="")
		{
			FeatureXMLFile().store(getStringOption_("unmatched_out"),unmatched_features);
		}
		
		//write pairs with too much deviation of the ratio
		if (getStringOption_("ratio_out")!="")
		{
			FeatureXMLFile().store(getStringOption_("ratio_out"),ratio_pairs);
		}
		
		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TOPPFFEVal tool;
	return tool.main(argc,argv);
}

/// @endcond




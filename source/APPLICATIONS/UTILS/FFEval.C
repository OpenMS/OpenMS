// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

using namespace OpenMS;
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
		registerDoubleOption_("rt_tol","<tol>",50.0,"Maximum allowed retention time deviation",false);
		registerDoubleOption_("mz_tol","<tol>",0.1,"Maximum allowed m/z deviation",false);

		addText_("Output options");
		registerOutputFile_("manual_out","<file>","","'manual' feature converted to featureXML",false);
		registerOutputFile_("intensity_out","<file>","","Gnuplot file of matched intensities",false);
		registerOutputFile_("ratio_out","<file>","","Gnuplot file of matched intensity ratios",false);
	}

	void medianError( vector<DoubleReal> a, vector<DoubleReal> b, DoubleReal& median, DoubleReal& q1, DoubleReal& q3)
	{
		vector<DoubleReal> errors;
		for (UInt i=0; i< a.size(); ++i) errors.push_back(fabs(b[i] - a[i]));
		sort(errors.begin(),errors.end());
		q1 = errors[errors.size()/4];
		median = errors[errors.size()/2];
		q3 = errors[(3*errors.size())/4];
	}

	DoubleReal meanError( vector<DoubleReal> a, vector<DoubleReal> b)
	{
		DoubleReal sum;
		for (UInt i=0; i< a.size(); ++i) sum += fabs(b[i] - a[i]);
		return  sum/a.size();
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
		
		//compare seek manual feature in automatic feature map
		DoubleReal mz_tol = getDoubleOption_("mz_tol");
		DoubleReal rt_tol = getDoubleOption_("rt_tol");
		UInt matched_single = 0;
		UInt matched_multi = 0;
		UInt matched_pairs = 0;
		vector<DoubleReal> a_int, m_int, a_ratio, m_ratio;
		Feature last_best_match;
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
				if ( f_a.getIntensity()>0.0 
					&& f_a.getRT()< (f_m.getRT() + rt_tol) 
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
						//cout << " - match: " << f_a.getRT() << "/" << f_a.getMZ() << " score: " <<  score << " -- best" << endl;				
					}
					else
					{
						//cout << " - match: " << f_a.getRT() << "/" << f_a.getMZ() << " score: " <<  score << endl;				
					}
				}
			}
			if (match_count==1) ++matched_single;
			if (match_count>1) ++matched_multi;
			
			//calculate statistics with the best match only
			if (match_count>=1)
			{
				a_int.push_back(best_match.getIntensity());
				m_int.push_back(f_m.getIntensity());
				
				//found pair
				if (Math::isOdd(m) && last_best_match!=Feature())
				{
					++matched_pairs;
					a_ratio.push_back(best_match.getIntensity()/last_best_match.getIntensity());
					m_ratio.push_back(f_m.getIntensity()/features_manual[m-1].getIntensity());
				}
			}
			
			//store last best match (for pairs)
			last_best_match = best_match;
		}
		
		cout << endl << endl;
		cout << "feature detection statistics: " << endl;
		cout << "  manual features: " << features_manual.size() << endl;
		cout << "  matches: " << matched_single + matched_multi << " (" << (100.0*(matched_single + matched_multi)/features_manual.size()) << "%)" << endl;
		cout << "  one match: " << matched_single << " (" << (100.0*matched_single/features_manual.size()) << "%)" << endl;
		cout << "  multiple matches: " << matched_multi << " (" << (100.0*matched_multi/features_manual.size()) << "%)" << endl;
		cout << "  intensity correlation: " << Math::pearsonCorrelationCoefficient(m_int.begin(),m_int.end(),a_int.begin(),a_int.end()) << endl;

		cout << endl << endl;
		cout << "pair detection statistics: " << endl;
		cout << "  manual pairs: " <<0.5*features_manual.size() << endl;
		cout << "  found: " << matched_pairs << " (" << (100.0*(matched_pairs)/(0.5*features_manual.size())) << "%)" << endl;
		DoubleReal median, q1, q3;
		medianError(m_ratio,a_ratio,median,q1,q3);
		cout << "  intensity ratio error q1/median/q3: " << q1 << " " << median << " " << q3 << endl;
		cout << "  intensity ratio error mean: " << meanError(m_ratio,a_ratio) << endl;
		
		
		//write intensity pair file
		if (getStringOption_("intensity_out")!="")
		{
			TextFile int_out;
			for (UInt i=0; i< a_int.size(); ++i)
			{
				int_out.push_back(String(m_int[i]) + " " + a_int[i]);		
			}
			int_out.store(getStringOption_("intensity_out"));
		}

		//write intensity pair file
		if (getStringOption_("ratio_out")!="")
		{
			TextFile ratio_out;
			for (UInt i=0; i< a_ratio.size(); ++i)
			{
				ratio_out.push_back(String(m_ratio[i]) + " " + a_ratio[i]);		
			}
			ratio_out.store(getStringOption_("ratio_out"));
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




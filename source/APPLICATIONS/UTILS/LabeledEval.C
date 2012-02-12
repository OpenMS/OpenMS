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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_LabeledEval LabeledEval
	
	@brief Evaluation tool for isotope-labeled quantitation experiments.

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_LabeledEval.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPLabeledEval
	: public TOPPBase
{
 public:
	TOPPLabeledEval()
		: TOPPBase("LabeledEval"," Evaluation tool for isotope-labeled quantitation experiments.",false)
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		addText_("Input options");
		registerInputFile_("in","<file>","","Feature result file");
		setValidFormats_("in", StringList::create("featureXML"));
		registerInputFile_("truth","<file>","","Expected result file.");
		setValidFormats_("truth", StringList::create("consensusXML"));
		registerDoubleOption_("rt_tol","<tol>",20.0,"Maximum allowed retention time deviation",false);
		registerDoubleOption_("mz_tol","<tol>",0.25,"Maximum allowed m/z deviation (divided by charge)",false);
	}

	String fiveNumbers(vector<DoubleReal> a, UInt decimal_places)
	{
		sort(a.begin(),a.end());
		return String::number(a[0],decimal_places) + " " + String::number(a[a.size()/4],decimal_places) + " " + String::number(a[a.size()/2],decimal_places) + " " + String::number(a[(3*a.size())/4],decimal_places) + " " + String::number(a.back(),decimal_places);
	}
	
	String fiveNumberQuotients(vector<DoubleReal> a, vector<DoubleReal> b, UInt decimal_places)
	{
		vector<DoubleReal> errors;
		for (Size i=0; i< a.size(); ++i) errors.push_back(a[i] / b[i]);
		return fiveNumbers(errors, decimal_places);
	}

	ExitCodes main_(int , const char**)
	{
		//load input features
		FeatureMap<> input;
		FeatureXMLFile().load(getStringOption_("in"),input);
		
		//load truth consensusXML
		ConsensusMap truth;
		ConsensusXMLFile().load(getStringOption_("truth"), truth);
		
		//parameters
		DoubleReal mz_tol = getDoubleOption_("mz_tol");
		DoubleReal rt_tol = getDoubleOption_("rt_tol");
		
		//seek manual feature in automatic feature map
		UInt matched_pairs = 0;
		UInt half_matched_pairs = 0;
		vector<DoubleReal> t_ratio, i_ratio, rt_diffs, mz_diffs;
		for (Size t=0; t<truth.size(); ++t)
		{
			if (truth[t].size()!=2)
			{
				cerr << "Error: consensus feature must contain exactly two elements!" << endl;
				continue;
			}
			vector<Feature> best_matches(2);
			vector<UInt> match_counts(2,0);
			vector<Peak2D> elements(2);
			elements[0] = *(truth[t].getFeatures().begin());
			elements[1] = *(++(truth[t].getFeatures().begin()));
			DoubleReal mz_tol_charged = mz_tol / truth[t].getCharge();
			for (Size e=0; e<2; ++e)
			{
				DoubleReal best_score = 0.0;
				for (Size i=0; i<input.size(); ++i)
				{
					const Feature& f_i = input[i];
					if ( fabs(f_i.getRT()-elements[e].getRT())<rt_tol
						&& fabs(f_i.getMZ()-elements[e].getMZ())<mz_tol_charged )
					{
						++match_counts[e];
						DoubleReal score = (1.0-fabs(f_i.getMZ()-elements[e].getMZ())/mz_tol_charged) * (1.0-fabs(f_i.getRT()-elements[e].getRT())/rt_tol);
						if (score > best_score)
						{
							best_score = score;
							best_matches[e] = f_i;
						}
					}
				}
			}

			//not matched
			if (match_counts[0]==0 && match_counts[1]==0)
			{
			}			
			//half matched
			else if ((match_counts[0]>0 && match_counts[1]==0) || (match_counts[0]==0 && match_counts[1]>0))
			{
				++half_matched_pairs;
			}
			//matched
			else
			{
				++matched_pairs;
				DoubleReal a_r = best_matches[0].getIntensity()/best_matches[1].getIntensity();
				t_ratio.push_back(a_r);
				DoubleReal m_r = elements[0].getIntensity()/elements[1].getIntensity();
				i_ratio.push_back(m_r);
				rt_diffs.push_back(best_matches[1].getRT()-best_matches[0].getRT());
				mz_diffs.push_back((best_matches[1].getMZ()-best_matches[0].getMZ())*truth[t].getCharge());
			}
		}
		
		cout << endl;
		cout << "pair detection statistics:" << endl;
		cout << "==========================" << endl;
		cout << "truth pairs: " << truth.size() << endl;
		cout << "input features: " << input.size() << endl;
		cout << endl;
		cout << "found: " << matched_pairs << " (" << String::number(100.0*matched_pairs/truth.size(),2) << "%)" << endl;
		cout << "half found : " << half_matched_pairs << " (" << String::number(100.0*half_matched_pairs/truth.size(),2) << "%)" << endl;
		cout << "not found : " << truth.size() - (matched_pairs+half_matched_pairs) << " (" << String::number(100.0-100.0*(matched_pairs+half_matched_pairs)/truth.size(),2) << "%)" << endl;
		cout << endl;
		cout << "relative pair ratios: " << fiveNumberQuotients(i_ratio,t_ratio,3) << endl;
		cout << "pair distance RT : " << fiveNumbers(rt_diffs, 2) << endl;
		cout << "pair distance m/z: " << fiveNumbers(mz_diffs, 2) << endl;
		
		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TOPPLabeledEval tool;
	return tool.main(argc,argv);
}

/// @endcond




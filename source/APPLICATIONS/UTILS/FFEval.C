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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FFEval FFEval
	
	@brief Evaluation tool for for feature detection algorithms.
		
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_FFEval.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFFEVal
	: public TOPPBase
{
 public:
	TOPPFFEVal()
		: TOPPBase("FFEVal","Evaluation tool for feature detection algorithms.",false)
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		addText_("Input options");
		registerInputFile_("in","<file>","","Feature input file, which contains the data to be tested against the truth file.");
		setValidFormats_("in", StringList::create("featureXML"));
		registerInputFile_("truth","<file>","","Truth feature file that defines what features should be found.");
		setValidFormats_("truth", StringList::create("featureXML"));
		registerOutputFile_("out","<file>","","Feature output file. If given, an annotated input file is written.", false);
		setValidFormats_("out", StringList::create("featureXML"));
		registerInputFile_("abort_reasons","<file>","","Feature file containing seeds with abort reasons.",false);
		setValidFormats_("abort_reasons", StringList::create("featureXML"));
		registerDoubleOption_("rt_tol","<double>",0.15,"Allowed tolerance of RT relative to average feature RT span.", false);
		setMinFloat_("rt_tol",0);
		registerDoubleOption_("rt_tol_abs","<double>",-1.0,"Allowed absolute tolerance of RT (overwrites 'rt_tol' if set above zero).", false);
		setMinFloat_("rt_tol_abs",-1);
		registerDoubleOption_("mz_tol","<double>",0.25,"Allowed tolerance in m/z (is divided by charge).", false);
		setMinFloat_("mz_tol",0);
	}
	
	/// Counts the number of features with meta value @p name equal to @p value
	UInt count(const FeatureMap<>& map, const String& name, const String& value="")
	{
		UInt count =0;
		for (Size i=0; i<map.size(); ++i)
		{
			if (map[i].metaValueExists(name))
			{
				if (value=="")
				{
					++count;
				}
				else
				{
					if (map[i].getMetaValue(name).toString()==value)
					{
						++count;
					}
				}
			} 
		}
		return count;
	}
	
	/// Returns the total number and percentage in parentheses
	String percentage(Size count, Size size)
	{
		return String(" (") + String::number(100.0*count/size,2) + "%)";
	}

	String fiveNumbers(vector<DoubleReal> a, UInt decimal_places)
	{
		sort(a.begin(),a.end());
		return String::number(a[0],decimal_places) + " " + String::number(a[a.size()/4],decimal_places) + " " + String::number(a[a.size()/2],decimal_places) + " " + String::number(a[(3*a.size())/4],decimal_places) + " " + String::number(a.back(),decimal_places);
	}

	ExitCodes main_(int , const char**)
	{
		//load data
		FeatureMap<> features_in,features_truth;
		FeatureXMLFile().load(getStringOption_("in"),features_in);
		features_in.sortByPosition();
		FeatureXMLFile().load(getStringOption_("truth"),features_truth);
		features_truth.sortByPosition();
		FeatureMap<> abort_reasons;
		if(getStringOption_("abort_reasons")!="")
		{
			FeatureXMLFile().load(getStringOption_("abort_reasons"),abort_reasons);
		}
		DoubleReal mz_tol = getDoubleOption_("mz_tol");
		writeDebug_(String("Final MZ tolerance: ") + mz_tol, 1);
		
		//determine average RT tolerance:
		//median feature RT span times given factor
		vector<DoubleReal> rt_spans;
		for (Size t=0; t<features_in.size(); ++t)
		{
			if (features_in[t].getConvexHulls().size()!=0)
			{
				rt_spans.push_back(features_in[t].getConvexHull().getBoundingBox().width());
			}
		}
		//feature convex hulls are available => relative RT span
		DoubleReal rt_tol = getDoubleOption_("rt_tol_abs");
		if (rt_tol<0.0)
		{
			if (rt_spans.size()!=0)
			{
				sort(rt_spans.begin(), rt_spans.end());
				rt_tol = getDoubleOption_("rt_tol")*rt_spans[rt_spans.size()/2];
			}
			else
			{
				writeLog_("Error: Input features do not have convex hulls. You have to set 'rt_tol_abs'!");
				return ILLEGAL_PARAMETERS;
			}
		}
		writeDebug_(String("Final RT tolerance: ") + rt_tol, 1);
		
		//general statistics
		std::vector<DoubleReal> ints_t;
		std::vector<DoubleReal> ints_i;
		std::vector<DoubleReal> ints_found;
		std::vector<DoubleReal> ints_missed;
		Map<String,UInt> abort_strings;

		for (Size m=0; m<features_truth.size(); ++m)
		{
			Feature& f_t =  features_truth[m];
			UInt match_count = 0;
			bool correct_charge = false;
			bool exact_centroid_match = false;
			Feature last_match;
			for (Size a=0; a<features_in.size(); ++a)
			{
				const Feature& f_i =  features_in[a];
				//RT match
				if (fabs(f_i.getRT()-f_t.getRT())<rt_tol)
				{
					DoubleReal charge_mz_tol = mz_tol / f_t.getCharge();
					//Exact m/z match
					if (fabs(f_i.getMZ()-f_t.getMZ())<charge_mz_tol)
					{
						++match_count;
						last_match = f_i;
						exact_centroid_match = true;
						if(f_i.getCharge()==f_t.getCharge()) correct_charge = true;
					}
					//Centroid is one trace off, but still contained in the convex hull
					else if (f_i.getConvexHull().encloses(f_t.getPosition())
									 &&
									 (
									  fabs(f_i.getMZ()+1.0/f_t.getCharge()-f_t.getMZ())<charge_mz_tol
									  ||
									  fabs(f_i.getMZ()-1.0/f_t.getCharge()-f_t.getMZ())<charge_mz_tol
									 )
						      )
					{
						++match_count;
						last_match = f_i;
						if(f_i.getCharge()==f_t.getCharge()) correct_charge = true;
					}
				}
			}
			
			f_t.setMetaValue("matches",match_count);
			if (match_count==1)
			{
				//flag matched feature with addition information
				if (correct_charge)
				{
					f_t.setMetaValue("correct_charge",String("true"));
					f_t.setMetaValue("intensity_ratio",last_match.getIntensity()/f_t.getIntensity());
				}
				else
				{
					f_t.setMetaValue("correct_charge",String("false"));
				}
				
				if (exact_centroid_match)
				{
					f_t.setMetaValue("exact_centroid_match",String("true"));
				}
				else
				{
					f_t.setMetaValue("exact_centroid_match",String("false"));
				}
			}
			//evaluation of correct features only
			if (match_count==1 && correct_charge)
			{
				ints_t.push_back(f_t.getIntensity());
				ints_i.push_back(last_match.getIntensity());
				ints_found.push_back(f_t.getIntensity());
			}
			else
			{
				ints_missed.push_back(f_t.getIntensity());
				
				//look up the abort reason of the nearest seed
				DoubleReal best_score_ab = 0;
				String reason = "";
				for (Size b=0; b<abort_reasons.size(); ++b)
				{
					const Feature& f_ab =  abort_reasons[b];
					if ( fabs(f_ab.getRT() - f_t.getRT()) <= rt_tol
						&& fabs(f_ab.getMZ() - f_t.getMZ()) <= mz_tol) 
					{
						DoubleReal score = (1.0-fabs(f_ab.getMZ()-f_t.getMZ())/mz_tol) * (1.0-fabs(f_ab.getRT()-f_t.getRT())/rt_tol);
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
		}

		//------------------------ general statistics ------------------------
		cout << endl;
		cout << "general information:" << endl;
		cout << "====================" << endl;
		cout << "input features: " << features_in.size() << endl;		
		cout << "truth features: " << features_truth.size() << endl;		

		//------------------------ matches ------------------------
		cout << endl;
		cout << "feature matching statistics:" << endl;
		cout << "============================" << endl;
		Size tmp = count(features_truth,"matches","0");
		cout << "no match: " << tmp << percentage(tmp,features_truth.size()) << endl;
		tmp = count(features_truth,"matches","1");
		cout << "one match: " << tmp << percentage(tmp,features_truth.size()) << endl;
		tmp = count(features_truth,"correct_charge","true");
		cout << " - correct charge: " << tmp << percentage(tmp,features_truth.size()) << endl;
		tmp = count(features_truth,"exact_centroid_match","true");
		cout << " - exact centroid match: " << tmp << percentage(tmp,features_truth.size()) << endl;
		tmp = features_truth.size() - count(features_truth,"matches","0") - count(features_truth,"matches","1");
		cout << "multiple matches: " << tmp << percentage(tmp,features_truth.size()) << endl;
		if (abort_reasons.size())
		{
			cout << "reasons for unmatched features:" << endl;
			for (Map<String,UInt>::iterator it=abort_strings.begin(); it!=abort_strings.end(); ++it)
			{
				cout << " - " << String(it->second).fillLeft(' ',4) << ": " << it->first << endl;
			}
		}
		//------------------------ intensity ------------------------
		cout << endl;
		cout << "intensity statistics:" << endl;
		cout << "=====================" << endl;
		if (ints_i.size()==0)
		{
			cout << "correlation of found features: nan" << endl;
			cout << "intensity distribution of found: 0.0 0.0 0.0 0.0 0.0" << endl;
			cout << "intensity distribution of missed: 0.0 0.0 0.0 0.0 0.0" << endl;
		}
		else
		{
			cout << "correlation of found features: " << pearsonCorrelationCoefficient(ints_i.begin(),ints_i.end(),ints_t.begin(),ints_t.end()) << endl;
			cout << "intensity distribution of found: " << fiveNumbers(ints_found,1) << endl;
			cout << "intensity distribution of missed: " << fiveNumbers(ints_missed,1) << endl;
		}
		
		//------------------------ charges ------------------------
		cout << endl;
		cout << "charge matches statistics:" << endl;
		cout << "===========================" << endl;
		Map<UInt,UInt> present_charges,found_charges;
		for (Size i=0;i<features_truth.size(); ++i)
		{
			UInt charge = features_truth[i].getCharge();
			present_charges[charge]++;
			if (features_truth[i].getMetaValue("correct_charge").toString()=="true")
			{
				found_charges[charge]++;
			}
		}
 		for (Map<UInt,UInt>::const_iterator it=present_charges.begin(); it!=present_charges.end(); ++it)
 		{
 			cout << "charge " << it->first << ": " << found_charges[it->first] << "/" << it->second << percentage(found_charges[it->first],it->second) << endl;
 		}
		
		//write output
		if (getStringOption_("out")!="")
		{
			FeatureXMLFile().store(getStringOption_("out"),features_truth);
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




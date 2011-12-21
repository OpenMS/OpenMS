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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_FFEval FFEval
	
	@brief Evaluation tool for feature detection algorithms.
		

  To plot the ROC curve you might use:

@code
  d = read.table("data.roc", skip=1, sep="\t")
  plot(d[,3],d[,4], xlim=c(0,1),ylim=c(0,1), xlab="FDR",ylab="TPR",main="ROC with varying intensity")
  lines(c(0,1),c(0,1))
@endcode

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_FFEval.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFFEval
	: public TOPPBase
{
 public:
	TOPPFFEval()
		: TOPPBase("FFEval","Evaluation tool for feature detection algorithms.",false)
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
		registerDoubleOption_("rt_tol","<double>",0.3,"Allowed tolerance of RT relative to average feature RT span.", false);
		setMinFloat_("rt_tol",0);
		registerDoubleOption_("rt_tol_abs","<double>",-1.0,"Allowed absolute tolerance of RT (overwrites 'rt_tol' if set above zero).", false);
		setMinFloat_("rt_tol_abs",-1);
		registerDoubleOption_("mz_tol","<double>",0.25,"Allowed tolerance in m/z (is divided by charge).", false);
		setMinFloat_("mz_tol",0);
		registerOutputFile_("out","<file>","","Feature output file. If given, an annotated input file is written.", false);
		setValidFormats_("out", StringList::create("featureXML"));
		registerInputFile_("abort_reasons","<file>","","Feature file containing seeds with abort reasons.",false);
		setValidFormats_("abort_reasons", StringList::create("featureXML"));
		registerOutputFile_("out_roc","<file>","","If given, a ROC curve file is created (ROC points based on intensity threshold)", false);
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
      if ( !rt_spans.empty() )
			{
				sort(rt_spans.begin(), rt_spans.end());
				rt_tol = getDoubleOption_("rt_tol")*rt_spans[rt_spans.size()/2];
			}
      else if (features_in.empty())
      {
        // do nothing, rt_tol does not really matter, as we will not find a match anyway, but we want to have the stats 
        // at the end, so we do not abort
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
			Size last_match_index = features_in.size()+1;
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
						exact_centroid_match = true;
						if(f_i.getCharge()==f_t.getCharge()) correct_charge = true;
						last_match_index = a;
					}
					//Centroid is one trace off, but still contained in the convex hull
					else if (f_i.getConvexHull().getBoundingBox().encloses(f_t.getPosition())
									 &&
									 (
									  fabs(f_i.getMZ()+1.0/f_t.getCharge()-f_t.getMZ())<charge_mz_tol
									  ||
									  fabs(f_i.getMZ()-1.0/f_t.getCharge()-f_t.getMZ())<charge_mz_tol
									 )
						      )
					{
						++match_count;
						last_match_index = a;
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
					f_t.setMetaValue("intensity_ratio",features_in[last_match_index].getIntensity()/f_t.getIntensity());
					features_in[last_match_index].setMetaValue("correct_hit", "true"); //flag the feature for ROC curve
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
				ints_i.push_back(features_in[last_match_index].getIntensity());
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
		Size no_match = count(features_truth,"matches","0");
		cout << "no match: " << no_match << percentage(no_match,features_truth.size()) << endl;
		Size one_match = count(features_truth,"matches","1");
		cout << "one match: " << one_match << percentage(one_match,features_truth.size()) << endl;
		Size charge_match = count(features_truth,"correct_charge","true");
		cout << " - correct charge: " << charge_match << percentage(charge_match,features_truth.size()) << endl;
		Size centroid_match = count(features_truth,"exact_centroid_match","true");
		cout << " - exact centroid match: " << centroid_match << percentage(centroid_match,features_truth.size()) << endl;
		Size multi_match = features_truth.size() - count(features_truth,"matches","0") - count(features_truth,"matches","1");
		cout << "multiple matches: " << multi_match << percentage(multi_match,features_truth.size()) << endl;
		Size incorrect_match = multi_match + one_match - charge_match;
		cout << "incorrect matches: " << incorrect_match << percentage(incorrect_match,features_truth.size()) << endl;
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
		if (ints_i.empty())
		{
			cout << "correlation of found features: nan" << endl;
		}
		else
		{
			cout << "correlation of found features: " << pearsonCorrelationCoefficient(ints_i.begin(),ints_i.end(),ints_t.begin(),ints_t.end()) << endl;
		}
		if (ints_found.empty())
		{
			cout << "intensity distribution of found: 0.0 0.0 0.0 0.0 0.0" << endl;
		}
		else
		{
			cout << "intensity distribution of found: " << fiveNumbers(ints_found,1) << endl;
		}
		if (ints_missed.empty())
		{
			cout << "intensity distribution of missed: 0.0 0.0 0.0 0.0 0.0" << endl;
		}
		else
		{
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
		
		//ROC curve
		if (getStringOption_("out_roc")!="")
		{
			TextFile tf;
			tf.push_back("false	correct	FDR	TPR");
			
			features_in.sortByIntensity(true);
			UInt f_correct = 0;
			UInt f_false = 0;
			DoubleReal found = features_in.size();
			DoubleReal correct = features_truth.size();
			for (Size i=0; i<features_in.size(); ++i)
			{
				if (features_in[i].metaValueExists("correct_hit"))
				{
					++f_correct;
				}
				else
				{
					++f_false;
				}
				tf.push_back(String(f_false) + "	" + f_correct + "	" + String::number(f_false/found,3) + "	" + String::number(f_correct/correct,3));
			}
			tf.store(getStringOption_("out_roc"));
		}
		
		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TOPPFFEval tool;
	return tool.main(argc,argv);
}

/// @endcond




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
	
	@brief Evaluation tool for isotope-labeled quantitation experiments.
	
	@todo Handling of truth convex hulls (Marc)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFFEVal
	: public TOPPBase
{
 public:
	TOPPFFEVal()
		: TOPPBase("FFEVal","Evaluation tool for isotope-labeled quantitation experiments.")
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
		registerDoubleOption_("rt_tol","<double>",0.15,"Allowed tolerance of RT relative to feature RT span.", false);
		setMinFloat_("rt_tol",0);
		setMaxFloat_("rt_tol",1);
		registerDoubleOption_("mz_tol","<double>",0.25,"Allowed tolerance in m/z (is devided by charge).", false);
		setMinFloat_("mz_tol",0);
		setMaxFloat_("mz_tol",1);
	}
	
	/// Counts the number of features with meta value @p name equal to @p value
	UInt count(const FeatureMap<>& map, const String& name, const String& value="")
	{
		UInt count =0;
		for (UInt i=0; i<map.size(); ++i)
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
	String percentage(UInt count, UInt size)
	{
		return String(" (") + String::number(100.0*count/size,2) + "%)";
	}
	
	ExitCodes main_(int , const char**)
	{
		//load data
		FeatureMap<> features_in,features_truth;
		FeatureXMLFile().load(getStringOption_("in"),features_in);
		features_in.sortByPosition();
		FeatureXMLFile().load(getStringOption_("truth"),features_truth);
		features_truth.sortByPosition();
		
		//general statistics
		std::vector<DoubleReal> ints_t;
		std::vector<DoubleReal> ints_i;
			
		for (UInt m=0; m<features_truth.size(); ++m)
		{
			Feature& f_t =  features_truth[m];
			UInt match_count = 0;
			bool correct_charge = false;
			bool exact_centroid_match = false;
			Feature last_match;
			for (UInt a=0; a<features_in.size(); ++a)
			{
				const Feature& f_i =  features_in[a];
				//RT match
				DoubleReal rt_tol = getDoubleOption_("rt_tol")*f_i.getConvexHull().getBoundingBox().width();
				if (fabs(f_i.getRT()-f_t.getRT())<rt_tol)
				{
					DoubleReal mz_tol = getDoubleOption_("mz_tol") / f_t.getCharge();
					//Exact m/z match
					if (fabs(f_i.getMZ()-f_t.getMZ())<mz_tol)
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
									  fabs(f_i.getMZ()+1.0/f_t.getCharge()-f_t.getMZ())<mz_tol
									  ||
									  fabs(f_i.getMZ()-1.0/f_t.getCharge()-f_t.getMZ())<mz_tol
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
					//intensity correlation
					ints_t.push_back(f_t.getIntensity());
					ints_i.push_back(last_match.getIntensity());
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
		UInt tmp = count(features_truth,"matches","0");
		cout << "no match: " << tmp << percentage(tmp,features_truth.size()) << endl;
		tmp = count(features_truth,"matches","1");
		cout << "one match: " << tmp << percentage(tmp,features_truth.size()) << endl;
		tmp = count(features_truth,"correct_charge","true");
		cout << " - correct charge: " << tmp << percentage(tmp,features_truth.size()) << endl;
		tmp = count(features_truth,"exact_centroid_match","true");
		cout << " - exact centroid match: " << tmp << percentage(tmp,features_truth.size()) << endl;
		tmp = features_truth.size() - count(features_truth,"matches","0") - count(features_truth,"matches","1");
		cout << "multiple matches: " << tmp << percentage(tmp,features_truth.size()) << endl;

		//------------------------ intensity ------------------------
		cout << endl;
		cout << "intensity statistics:" << endl;
		cout << "=====================" << endl;
		cout << "correlation of correct features: " << pearsonCorrelationCoefficient(ints_i.begin(),ints_i.end(),ints_t.begin(),ints_t.end()) << endl;
		
		//------------------------ charges ------------------------
		cout << endl;
		cout << "charge matches statistics:" << endl;
		cout << "===========================" << endl;
		Map<UInt,UInt> present_charges,found_charges;
		for (UInt i=0;i<features_truth.size(); ++i)
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




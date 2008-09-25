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
		for (UInt m=0; m<features_truth.size(); ++m)
		{
			Feature& f_t =  features_truth[m];
			UInt match_count = 0;
			bool charge_ok = false;
			Feature best_match;
			for (UInt a=0; a<features_in.size(); ++a)
			{
				const Feature& f_i =  features_in[a];
				if (f_i.getConvexHull().encloses(f_t.getPosition()))
				{
					++match_count;
					if(f_i.getCharge()==f_t.getCharge()) charge_ok = true;
				}
				else 
				{
					if (fabs(f_i.getRT()-f_t.getRT())<0.25*f_i.getConvexHull().getBoundingBox().width())
					{
						if (f_t.getCharge()!=0 && fabs(f_i.getMZ()-f_t.getMZ())<0.25/f_t.getCharge())
						{
							++match_count;
							if(f_i.getCharge()==f_t.getCharge()) charge_ok = true;
						}
					}
				}
			}
			
			if (match_count==0)
			{
				f_t.setMetaValue("matches",String("none"));
			}
			else if (match_count==1)
			{
				f_t.setMetaValue("matches",String("one"));
				if (charge_ok)
				{
					f_t.setMetaValue("charge_correct",String("true"));
				}
				else
				{
					f_t.setMetaValue("charge_correct",String("false"));
				}
			}
			else if (match_count>1)
			{
				f_t.setMetaValue("matches",String("mutiple"));
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
		UInt tmp = count(features_truth,"matches","none");
		cout << "no match: " << tmp << percentage(tmp,features_truth.size()) << endl;
		tmp = count(features_truth,"matches","one");
		cout << "one match: " << tmp << percentage(tmp,features_truth.size()) << endl;
		tmp = count(features_truth,"matches","multiple");
		cout << "multiple matches: " << tmp << percentage(tmp,features_truth.size()) << endl;
		
		//------------------------ charges ------------------------
		cout << endl;
		cout << "charge matches statistics:" << endl;
		cout << "===========================" << endl;
		tmp = count(features_truth,"charge_correct","true");
		cout << "one match and correct charge: " << tmp << percentage(tmp,features_truth.size()) << endl;
		Map<UInt,UInt> present_charges,found_charges;
		for (UInt i=0;i<features_truth.size(); ++i)
		{
			UInt charge = features_truth[i].getCharge();
			present_charges[charge]++;
			if (features_truth[i].getMetaValue("charge_correct").toString()=="true")
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




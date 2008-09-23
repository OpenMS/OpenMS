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
		: TOPPBase("FFEVal"," Evaluation tool for isotope-labeled quantitation experiments.")
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		addText_("Input options");
		registerInputFile_("features","<file>","","Feature result file");
		setValidFormats_("features", StringList::create("featureXML"));
		registerInputFile_("manual","<file>","","Manual result file");
		setValidFormats_("manual", StringList::create("featureXML"));
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
		//load data
		FeatureMap<> features_automatic,features_manual;
		FeatureXMLFile().load(getStringOption_("features"),features_automatic);
		features_automatic.sortByPosition();
		FeatureXMLFile().load(getStringOption_("manual"),features_manual);
		features_manual.sortByPosition();
		
		//general statistics
		UInt matched_single=0;
		UInt matched_multi=0;
				
		for (UInt m=0; m<features_manual.size(); ++m)
		{
			const Feature& f_m =  features_manual[m];
			UInt match_count = 0;
			Feature best_match;
			for (UInt a=0; a<features_automatic.size(); ++a)
			{
				const Feature& f_a =  features_automatic[a];
				if (f_a.getConvexHull().encloses(f_m.getPosition()))
				{
					++match_count;
				}
			}
			if (match_count==1) ++matched_single;
			if (match_count>1) ++matched_multi;
		}
		
		cout << endl;
		cout << "feature detection statistics:" << endl;
		cout << "=============================" << endl;
		cout << "  manual features: " << features_manual.size() << endl;
		cout << "  matches: " << matched_single + matched_multi << " (" << String::number(100.0*(matched_single + matched_multi)/features_manual.size(),2) << "%)" << endl;
		cout << "    one match: " << matched_single << " (" << String::number(100.0*matched_single/features_manual.size(),2) << "%)" << endl;
		cout << "    multiple matches: " << matched_multi << " (" << String::number(100.0*matched_multi/features_manual.size(),2) << "%)" << endl;
		
		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TOPPFFEVal tool;
	return tool.main(argc,argv);
}

/// @endcond




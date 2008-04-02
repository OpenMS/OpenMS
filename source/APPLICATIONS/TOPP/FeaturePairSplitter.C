// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeaturePairsXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <map>
#include <iostream>
#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------


/**
	@page FeaturePairSplitter FeaturePairSplitter

	@brief Splits a featurePairsXML file into two featureXML files.

	This is just a small utility.  The features are copied from the pairs.  The
	relative order of features is preserved.  For example, the first two features
	of each output file belong to each other, then the second two, and so on.
	The <i>quality</i> information of the feature pairs can be written to a third
	file.

  A typical file name extension for the input would be ".featurePairsXML".

 	A typical file name extension for the two output files would be '.featureXML'.

	The qualities are written one per line; a typical file name extension would be '.txt'.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeaturePairSplitter
      : public TOPPBase
{
 public:
  TOPPFeaturePairSplitter()
		: TOPPBase("FeaturePairSplitter","split a feature pairs file into two featureXML files and a qualities file")
  {
  }

 protected:
  typedef std::vector < ElementPair < Feature > > FeaturePairVector;
  
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in","<file>","","feature pair file ");
		setValidFormats_("in",StringList::create("FeaturePairsXML"));
    registerOutputFile_("out1","<file>","","first FeatureXML output file ",false);
	  setValidFormats_("out1",StringList::create("FeatureXML"));
    registerOutputFile_("out2","<file>","","second FeatureXML output file ",false);
	  setValidFormats_("out2",StringList::create("FeatureXML"));
    registerOutputFile_("qual","<file>","","pair qualtities output file",false);
  }

  ExitCodes main_(int , const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

		// file names
		String in = getStringOption_("in");
		String out1 = getStringOption_("out1");
		bool const write_out1 = !out1.empty();
		String out2 = getStringOption_("out2");
		bool const write_out2 = !out2.empty();
		String qual = getStringOption_("qual");
		bool const write_qual = !qual.empty();

		// load data from input file.
    FeaturePairVector feature_pairs;
		FeaturePairsXMLFile feature_pairs_file;
		feature_pairs_file.load(in,feature_pairs);

		// store the data
		FeatureMap<> first_feature_map, second_feature_map;
		vector<double> qualities_vector;
		for ( FeaturePairVector::const_iterator iter = feature_pairs.begin();
					iter != feature_pairs.end();
					++iter
				)
		{
			if ( write_out1 ) first_feature_map.push_back(iter->getFirst());
			if ( write_out2 ) second_feature_map.push_back(iter->getSecond());
			if ( write_qual ) qualities_vector.push_back(iter->getQuality());
		}

		// write the data to files
		FeatureXMLFile f;
		if ( write_out1 )
		{
			FeatureXMLFile f;
			f.store(out1,first_feature_map);
		}
		if ( write_out2 )
		{
			FeatureXMLFile f;
			f.store(out2,second_feature_map);
		}
		if ( write_qual )
		{
			ofstream qualities_file(qual.c_str());
			copy(qualities_vector.begin(),qualities_vector.end(), ostream_iterator<double>(qualities_file,"\n") );
		}

		return EXECUTION_OK;

	} // main_()

}; // TOPPFeaturePairSplitter


int main( int argc, const char** argv )
{
  TOPPFeaturePairSplitter tool;
  return tool.main(argc,argv);
}

/// @endcond

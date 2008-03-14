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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

#include<vector>
#include<algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TextExporter TextExporter
	
	@brief This application converts several OpenMS XML formats
	(namely featureXML, consensusXML and idXML) to text files.
	These text files can be easily read using other applications
	such as R/Matlab/Excel etc.
  
  @todo Add support for IdXML format (Andreas)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


typedef ConsensusFeature< > ConsensusFeatureType;

// Sorts consensus elements by size
struct ConsensusElementComparator
{

		inline bool operator() (const ConsensusFeatureType & x, const ConsensusFeatureType & y)
		{
			return x.size() < y.size();
		}
};


class TOPPTextExporter
	: public TOPPBase
{
	public:
		TOPPTextExporter()
			: TOPPBase("TextExporter","Exports various XML formats to a text file")
		{
			
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
      registerInputFile_("in","<file>","","input file");
      registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
			setValidStrings_("in_type",StringList::create("featureXML,consensusXML,idXML"));

      registerOutputFile_("out","<file>","","text file");
		}
	
		ExitCodes main_(int , const char**)
		{

			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			String in = getStringOption_("in");
			String out = getStringOption_("out");
        
      //input file type
      FileHandler fh;
      FileHandler::Type in_type = fh.nameToType(getStringOption_("in_type"));

      if (in_type==FileHandler::UNKNOWN)
      {
        in_type = fh.getTypeByFileName(in);
        writeDebug_(String("Input file type (from file extention): ") + fh.typeToName(in_type), 2);
      }
  
      if (in_type==FileHandler::UNKNOWN)
      {
        in_type = fh.getTypeByContent(in);
        writeDebug_(String("Input file type (from content): ") + fh.typeToName(in_type), 2);
      }
  
      if (in_type==FileHandler::UNKNOWN)
      {
        writeLog_("Error: Could not determine input file type!");
        return PARSE_ERROR;
      }
      
      if (in_type == FileHandler::FEATURE)
      {
  			 //-------------------------------------------------------------
        // loading input
        //-------------------------------------------------------------

        typedef FeatureMap<> FeatureMapType;
        FeatureMapType feature_map;
        FeatureXMLFile f;
        f.load(in,feature_map);             
  		
				 // text output
        ofstream outstr( out.c_str() );

				// stores one feature per line
				outstr << "# rt, mz, intensity, charge, overall_quality, rt_quality, mz_quality, rt_start, rt_end" << endl;
        for (FeatureMap< >::const_iterator citer = feature_map.begin();
             	citer != feature_map.end();
             ++citer)
        {
                outstr << citer->getPosition()[0] << " " << citer->getPosition()[1] << " " << citer->getIntensity();
                outstr << " " << citer->getCharge();
                outstr << " " << citer->getOverallQuality();
                outstr << " " << citer->getQuality(0) << " " << citer->getQuality(1);
                outstr << " " << citer->getConvexHulls().begin()->getBoundingBox().minX();
                outstr << " " << citer->getConvexHulls().begin()->getBoundingBox().maxX();
                outstr << endl;
        }
        outstr.close();
			
      }
      else if (in_type == FileHandler::CONSENSUSXML)
      {
					ConsensusMap< > cmap;
					vector<FeatureMap<> > feat_maps(100);
						
					// Yes, I know that this is ugly.  This is a problem with
					// the ConsensusMap, so don't bug me but speak to Eva. :-) 
					cmap.getMapVector().resize(100);
					for (UInt i = 0; i < 100; ++i)
					{
						cmap.getMapVector()[i] = &(feat_maps[i]);
					}

					/// No progress logging implemented for ConsensusXMLFile
					ConsensusXMLFile().load(in,cmap);
													
					UInt nr_conds = cmap.getFilenames().size();
					
					// A consensus feature map consisting of many feature maps will often
					// contain a lot of singleton features (i.e. features detected only in one
					// LC-MS map). We want to put these features at the end of the text file.
					// => sort consensus elements by size 
					sort(cmap.begin(),cmap.end(),ConsensusElementComparator() );
      
					ofstream txt_out( out.c_str() );
					
					// write header
					txt_out << "# consensus_rt consensus_mz ";
					for (UInt i=0;i<nr_conds;++i)
					{
						txt_out << "exp_" + String(i+1) + " ";
					}
					txt_out << endl;
					
					for (ConsensusMap< >::iterator cmap_it = cmap.begin(); cmap_it != cmap.end();++cmap_it)
					{
						// write consensus rt and m/z
						txt_out << cmap_it->getPosition()[0] << " " << cmap_it->getPosition()[1] << " ";
							
						UInt curr_cond = 0;						 		 																
						for ( ConsensusFeature< >::Group::const_iterator group_it = cmap_it->begin(); group_it != cmap_it->end(); ++group_it)
						{			
							UInt this_cond = group_it->getMapIndex();
							
							// print 0 (not available) for missing values
							while ( (curr_cond) != this_cond )
							{
								txt_out << "0 ";
								++curr_cond;
							}
						 		 																		    																																						
							txt_out << group_it->getElement().getIntensity() << " "; 
							++curr_cond;						
							
						}	// end for all features in this consensus element	    																																																					    			
						
						// append zeros for missing feature maps / conditions ( we start counting at zero)
						while (curr_cond <= ( nr_conds - 1) )
						{
							txt_out << "0 ";
							++curr_cond;
						}
						
						txt_out << endl;
				} // end for all elements in consensus map 
				
				txt_out.close();
      }
      else
      {
        writeLog_("Unknown input file type given. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;          
      }
			
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPTextExporter t;
	return t.main(argc,argv);
}

/// @endcond

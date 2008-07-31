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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

#include <vector>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TextExporter TextExporter
	
	@brief This application converts several %OpenMS XML formats
	(namely featureXML, consensusXML and idXML) to text files.
	These text files can be easily read using other applications
	such as R, Matlab, Excel, etc.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


// Sorts consensus elements by size
struct ConsensusElementComparator
{

		inline bool operator() (const ConsensusFeature & x, const ConsensusFeature & y)
		{
			return x.size() < y.size();
		}
};


class TOPPTextExporter
	: public TOPPBase
{
	public:
		TOPPTextExporter()
			: TOPPBase("TextExporter","Exports various XML formats to a text file.")
		{
			
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
      registerInputFile_("in","<file>","","input file ");
    	setValidFormats_("in",StringList::create("featureXML,consensusXML,idXML"));
			registerFlag_("proteins_only", "set this flag if you want only protein information from an idXML file");
			registerFlag_("peptides_only", "set this flag if you want only peptide information from an idXML file");

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
      FileHandler::Type in_type = FileHandler::getType(in);
      writeDebug_(String("Input file type: ") + FileHandler::typeToName(in_type), 2);
  
      if (in_type==FileHandler::UNKNOWN)
      {
        writeLog_("Error: Could not determine input file type!");
        return PARSE_ERROR;
      }
      
      if (in_type == FileHandler::FEATUREXML)
      {
  			 //-------------------------------------------------------------
        // loading input
        //-------------------------------------------------------------

        FeatureMap<> feature_map;
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

					if (citer->getConvexHulls().size() > 0)
					{
          	outstr << " " << citer->getConvexHulls().begin()->getBoundingBox().minX();
          	outstr << " " << citer->getConvexHulls().begin()->getBoundingBox().maxX();
					}
					else
					{
						outstr << " -1";
						outstr << " -1";
					}
          outstr << endl;
        }
        outstr.close();
			
      }
      else if (in_type == FileHandler::CONSENSUSXML)
      {
				ConsensusMap cmap;
				vector<FeatureMap<> > feat_maps(100);

				/// No progress logging implemented for ConsensusXMLFile
				ConsensusXMLFile().load(in,cmap);
												
				// A consensus feature map consisting of many feature maps will often
				// contain a lot of singleton features (i.e. features detected only in one
				// LC-MS map). We want to put these features at the end of the text file.
				// => sort consensus elements by size 
				sort(cmap.begin(),cmap.end(),ConsensusElementComparator() );
    
				ofstream txt_out( out.c_str() );
				
				//Write file descriptions
		 		const ConsensusMap::FileDescriptions& descs = cmap.getFileDescriptions();
	 			txt_out << "#Source file descriptions:" << endl;
	 			txt_out << "#" << endl;
		 		for (ConsensusMap::FileDescriptions::const_iterator it=descs.begin(); it!=descs.end(); ++it)
		 		{
					txt_out << "# identifier " << it->first << ": " << endl;
					txt_out << "#   filename : " << it->second.filename << endl;
					String label = it->second.label;
					label.trim();
					if (label!="") txt_out << "#   label : " << it->second.label << endl;
					if (it->second.size!=0) txt_out  << "#   size : " << it->second.size << endl;
					txt_out << "#" << endl;
		 		}
				
				// write header
				txt_out << "#consensus_rt	consensus_mz	consensus_intensity	quality	";
		 		for (ConsensusMap::FileDescriptions::const_iterator it=descs.begin(); it!=descs.end(); ++it)
		 		{
					txt_out << "	intensity_" << it->first;
				}
				txt_out << endl;
				
				for (ConsensusMap::iterator cmap_it = cmap.begin(); cmap_it != cmap.end();++cmap_it)
				{
					// write consensus rt and m/z
					txt_out << cmap_it->getPosition()[0] << "	" << cmap_it->getPosition()[1] << "	" << cmap_it->getIntensity() << "	" << cmap_it->getQuality();
					
					//determine present values	
					Map<UInt,DoubleReal> intensities;			 		 																
					for ( ConsensusFeature::HandleSetType::const_iterator group_it = cmap_it->begin(); group_it != cmap_it->end(); ++group_it)
					{
						intensities[group_it->getMapIndex()] = group_it->getIntensity();
					}
					
					//print all values (0.0 for missing ones)
					for (ConsensusMap::FileDescriptions::const_iterator it=descs.begin(); it!=descs.end(); ++it)
			 		{
						if (intensities.has(it->first))
						{
							txt_out << "	" << intensities[it->first];
						}
						else
						{
							txt_out << "	0.0";
						}
					}
					txt_out << endl;
				}
				txt_out.close();
      }
			else if (in_type == FileHandler::IDXML)
			{
				vector<ProteinIdentification> prot_ids;
				vector<PeptideIdentification> pep_ids;
				IdXMLFile().load(in, prot_ids, pep_ids);
				
				
				ofstream txt_out(out.c_str());

				for (vector<ProteinIdentification>::const_iterator it = prot_ids.begin(); it != prot_ids.end(); ++it)
				{
					String actual_id = it->getIdentifier();
					if (!getFlag_("peptides_only"))
					{
						// protein id header 
						txt_out << "# Run ID, Score Type, Score Direction, Date/Time, Search Engine Version " << endl;
						txt_out << actual_id << " " 
										<< it->getScoreType() << " ";
						if (it->isHigherScoreBetter())
						{
							txt_out << "higher-score-better ";
						}
						else
						{
							txt_out << "lower-score-better ";
						}
						// using ISODate ensures that TOPP tests will run through regardless of locale setting
						txt_out << it->getDateTime().toString(Qt::ISODate).toStdString() << " "
										<< it->getSearchEngineVersion() << endl;

						// search parameters
						ProteinIdentification::SearchParameters sp = it->getSearchParameters();
						txt_out << "# Search parameters of ID=" << actual_id 
										<< ": db=" << sp.db 
										<< ", db_version=" << sp.db_version 
										<< ", taxonomy=" << sp.taxonomy 
										<< ", charges=" << sp.charges 
										<< ", mass_type=";
						if (sp.mass_type == ProteinIdentification::MONOISOTOPIC)
						{
							txt_out << "monoisotopic";
						}
						else
						{
							txt_out << "average";
						}
						txt_out << ", fixed_modifications=";
						for (vector<String>::const_iterator mit = sp.fixed_modifications.begin(); mit != sp.fixed_modifications.end(); ++mit)
						{
							if (mit != sp.fixed_modifications.begin())
							{
								txt_out << ";";
							}
							txt_out << *mit;
						}
						txt_out << ", variable_modifications=";
						for (vector<String>::const_iterator mit = sp.variable_modifications.begin(); mit != sp.variable_modifications.end(); ++mit)
						{
							if (mit != sp.variable_modifications.begin())
							{
								txt_out << ";";
							}
							txt_out << *mit;
						}
						txt_out << ", enzyme=";
						switch (sp.enzyme)
						{
							case ProteinIdentification::TRYPSIN:
								txt_out << "Trypsin";
								break;
							case ProteinIdentification::PEPSIN_A:
								txt_out << "PepsinA";
								break;
							case ProteinIdentification::PROTEASE_K:
								txt_out << "ProteaseK";
								break;
							case ProteinIdentification::CHYMOTRYPSIN:
								txt_out << "ChymoTrypsin";
								break;
							default:
								txt_out << "unknown";
						}
						txt_out << ", missed_cleavages=" << sp.missed_cleavages 
										<< ", peak_mass_tolerance=" << sp.peak_mass_tolerance 
										<< ", precursor_mass_tolerance=" << sp.precursor_tolerance << endl;
						
								
						// header of protein hits
						txt_out << "# Protein Hits: Score, Rank, Accession, Sequence" << endl;
					
						for (vector<ProteinHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
						{
							txt_out << pit->getScore() << " " 
											<< pit->getRank() << " " 
											<< pit->getAccession() << " " 
											<< pit->getSequence() << endl;
						}
					}

					
					if (!getFlag_("proteins_only"))
					{
						// slight improvement on big idXML files with many different runs:
						// index the identifiers and peptide ids to avoid running over them 
						// again and again
						for (vector<PeptideIdentification>::const_iterator pit = pep_ids.begin(); pit != pep_ids.end(); ++pit)
						{
							if (pit->getIdentifier() == actual_id)
							{
								// header of peptide idenfication
								txt_out << "# RunID, RT, m/z, ScoreType, Score Direction" << endl;
								txt_out << actual_id << " ";
								
								if (pit->metaValueExists("RT"))
								{
									txt_out << (double)pit->getMetaValue("RT") << " ";
								}
								else
								{
									txt_out << "-1 ";
								}

								if (pit->metaValueExists("MZ"))
								{
									txt_out << (double)pit->getMetaValue("MZ") << " ";
								}
								else
								{
									txt_out << "-1 ";
								}
								
								txt_out	<< pit->getScoreType() << " ";
								if (pit->isHigherScoreBetter())
            		{
              		txt_out << "higher-score-better ";
            		}
            		else
            		{
              		txt_out << "lower-score-better ";
           			}
								txt_out << endl;

											
								// header of peptide hits
            		txt_out << "# Peptide Hits: Score, Rank, Sequence, Charge, AABefore, AAAfter, Accessions" << endl;

            		for (vector<PeptideHit>::const_iterator ppit = pit->getHits().begin(); ppit != pit->getHits().end(); ++ppit)
            		{
              		txt_out << ppit->getScore() << " "
                		      << ppit->getRank() << " "
                  		    << ppit->getSequence() << " "
													<< ppit->getCharge() << " "
													<< ppit->getAABefore() << " " 
													<< ppit->getAAAfter() << " ";
									for (vector<String>::const_iterator ait = ppit->getProteinAccessions().begin(); ait != ppit->getProteinAccessions().end(); ++ait)
									{
										if (ait != ppit->getProteinAccessions().begin())
										{
											txt_out << ";";
										}
										txt_out << *ait;
									}
									txt_out << endl;
								}
            	}
						}
					}
				}
				
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

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
	@page TOPP_TextExporter TextExporter
	
	@brief This application converts several %OpenMS XML formats
	(namely featureXML, consensusXML and idXML) to text files.
	These text files can be easily read using other applications
	such as R, Matlab, Excel, etc.
	
	@todo Add identifications output to featureXML and consensusXML and consider no_ids flag (Andreas, Clemens, Chris, Nico)

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_TextExporter.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{

	/// Wrapper class to implement printing of FeatureHandle
	struct FeatureHandlePrinter
	{
		FeatureHandlePrinter(const FeatureHandle &rhs) : ref_(rhs) {}
	 	const FeatureHandle &ref_;
	};

	// There is no standard how to print a nan (not-a-number) value.  So we do
	// this on our own.
	namespace 
	{
		const char nan[] = "nan"; // that's what Linux GCC uses, and gnuplot understands.

		template < typename NumberT >
		std::ostream & printValueOrNan( std::ostream & os, NumberT thing )
		{
			if ( !isnan(thing) )
			{
				return os << thing;
			}
			else
			{
				return os << nan;
			}
		}
	}

	/// Output operator for a FeatureHandlePrinter.
	std::ostream & operator << ( std::ostream& os, const FeatureHandlePrinter& rhs)
	{
		const UInt exponent_extra_digits = 6;
		const UInt charge_digits = 5;
		const unsigned prec_save = os.precision();
		os << std::setprecision(writtenDigits<FeatureHandle::CoordinateType>())
			 << std::setw(writtenDigits<FeatureHandle::CoordinateType>()+exponent_extra_digits);
		printValueOrNan(os,rhs.ref_.getRT());
		os << ' '
			 << std::setw(writtenDigits<FeatureHandle::CoordinateType>()+exponent_extra_digits);
		printValueOrNan(os,rhs.ref_.getMZ());
		os << ' '
			 << std::setprecision(writtenDigits<FeatureHandle::IntensityType>())
			 << std::setw(writtenDigits<FeatureHandle::IntensityType>()+exponent_extra_digits);
		printValueOrNan(os,rhs.ref_.getIntensity());
		os << ' '
			 << std::setw(charge_digits)
			 << rhs.ref_.getCharge()
			 << std::setprecision(prec_save);
		return os;
	}

	/// Wrapper class to implement printing of ConsensusFeature
	struct ConsensusFeaturePrinter
	{
		ConsensusFeaturePrinter(const ConsensusFeature &rhs) : ref_(rhs) {}
	 	const ConsensusFeature &ref_;
	};

	/// Output operator for a ConsensusFeaturePrinter.
	std::ostream & operator << ( std::ostream& os, const ConsensusFeaturePrinter& rhs)
	{
		const UInt exponent_extra_digits = 6;
		const UInt charge_digits = 5;
		const unsigned prec_save = os.precision();
		os << std::setprecision(writtenDigits<FeatureHandle::CoordinateType>())
			 << std::setw(writtenDigits<FeatureHandle::CoordinateType>()+exponent_extra_digits);
		printValueOrNan(os,rhs.ref_.getRT());
		os << ' '
			 << std::setw(writtenDigits<FeatureHandle::CoordinateType>()+exponent_extra_digits);
		printValueOrNan(os,rhs.ref_.getMZ());
		os << ' '
			 << std::setprecision(writtenDigits<FeatureHandle::IntensityType>())
			 << std::setw(writtenDigits<FeatureHandle::IntensityType>()+exponent_extra_digits);
		printValueOrNan(os,rhs.ref_.getIntensity());
		os << ' '
			 << std::setw(charge_digits)
			 << rhs.ref_.getCharge()
			 << std::setprecision(prec_save);
		return os;
	}

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
      registerInputFile_("in","<file>","","Input file ");
    	setValidFormats_("in",StringList::create("featureXML,consensusXML,idXML"));
      registerOutputFile_("out","<file>","","Output file. Only used for FeatureXML and IdXML.",false);
      registerFlag_("no_ids","Suppresses output of identification data for consensusXML and featureXML");
			addEmptyLine_();
			addText_("Options for IdXML files:");
			registerFlag_("proteins_only", "Set this flag if you want only protein information from an idXML file");
			registerFlag_("peptides_only", "Set this flag if you want only peptide information from an idXML file");
			registerFlag_("peptides_only_csv","Set this flag if you want only peptide information from an idXML file in csv format",false);
			addEmptyLine_();
			addText_("Options for ConsensusXML files:");
			registerOutputFile_("consensus_centroids","<file>","","Centroids of consensus features",false);
			registerOutputFile_("consensus_elements","<file>","","Elements of consensus features",false);
			registerOutputFile_("consensus_features","<file>","","Consensus features and contained elements from all maps (writes 'nan's if element is missing)",false);
			addText_("Each of the consensus_... files is created as requested.");
			registerStringOption_("sorting_method","<method>","none","Sorting method",false);
			setValidStrings_("sorting_method",StringList::create("none,RT,MZ,RT_then_MZ,intensity,quality_decreasing,quality_increasing"));
			registerFlag_("sort_by_maps","Apply a stable sort by the covered maps, lexicographically",false);
			registerFlag_("sort_by_size","Apply a stable sort by decreasing size (i.e., the number of elements)",false);
			addText_("Sorting options can be combined.  The precedence is: sort_by_size, sort_by_maps, sorting_method");
			registerFlag_("first_dim_rt","If this flag is set the first_dim RT of the peptide hits will also be printed (if present).");
			return;
		}
	
		/*

		An example for ConsensusXML:
	
		splot 'consensus_features.wsv' using 1:1:1 w l lw 2, 'consensus_features.wsv' using 1:1:1:($5-$1):(0):(0) with vectors nohead, 'consensus_features.wsv' using 1:1:1:(0):($9-$1):(0) with vectors nohead, 'consensus_features.wsv' using 1:1:1:(0):(0):($13-$1) with vectors nohead, 'consensus_features.wsv' using 5:1:1 w l, 'consensus_features.wsv' using 1:9:1 w l,'consensus_features.wsv' using 1:1:13 w l
	
		*/


		ExitCodes main_(int , const char**)
		{

			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			String in = getStringOption_("in");
			String out = getStringOption_("out");
			UInt counter = 0;
			bool without_header_repetition = getFlag_("peptides_only_csv");
      bool no_ids = getFlag_("no_ids");
      bool first_dim_rt = getFlag_("first_dim_rt");
        
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
				if ( out != "" )
				{
					writeLog_("Option 'out' is not functional for Consensusxml.  Use the 'consensus_...' options instead.");
					printUsage_();
					return ILLEGAL_PARAMETERS;          
				}

				String consensus_centroids = getStringOption_("consensus_centroids");
				String consensus_elements = getStringOption_("consensus_elements");
				String consensus_features = getStringOption_("consensus_features");
				String sorting_method = getStringOption_("sorting_method");
				bool sort_by_maps = getFlag_("sort_by_maps");
				bool sort_by_size = getFlag_("sort_by_size");

				ConsensusMap consensus_map;
				ConsensusXMLFile consensus_xml_file;

				consensus_xml_file.load(in,consensus_map);

				if ( sorting_method == "none" )
				{
					// don't sort in this case
				}
				else if ( sorting_method == "RT")
				{
					consensus_map.sortByRT();
				}
				else if ( sorting_method == "MZ")
				{
					consensus_map.sortByMZ();
				}
				else if ( sorting_method == "RT_then_MZ")
				{
					consensus_map.sortByPosition();
				}
				else if ( sorting_method == "intensity")
				{
					consensus_map.sortByIntensity();
				}
				else if ( sorting_method == "quality_decreasing")
				{
					consensus_map.sortByQuality(true);
				}
				else if ( sorting_method == "quality_increasing")
				{
					consensus_map.sortByQuality(false);
				}
				else
				{
					// wrong sorting methods should already be trapped by setValidStrings_, but who knows...
					writeLog_(String("Error: unknown sorting method: ")+sorting_method);
					return PARSE_ERROR;
				}
	
				if ( sort_by_maps )
				{
					consensus_map.sortByMaps();
				}

				if ( sort_by_size )
				{
					consensus_map.sortBySize();
				}

				DateTime date_time;
				date_time.now();
				String date_time_now = date_time.get();

				// ----------------------------------------------------------------------

				if ( !consensus_centroids.empty() )
				{
					std::ofstream consensus_centroids_file(consensus_centroids.c_str());
					if (!consensus_centroids_file)
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, consensus_centroids);
					}

					consensus_centroids_file
						<< "#  Centroids of consensus features extracted from " << in 
						<< " on " << date_time_now << std::endl
						<< "# RT MZ Intensity Charge" << std::endl;
					for ( ConsensusMap::const_iterator cmit = consensus_map.begin(); cmit != consensus_map.end(); ++cmit )
					{
						consensus_centroids_file << ConsensusFeaturePrinter(*cmit) << std::endl;
					}
					consensus_centroids_file.close();
				}

				// ----------------------------------------------------------------------
			
				if ( !consensus_elements.empty() )
				{
					std::ofstream consensus_elements_file(consensus_elements.c_str());
					if (!consensus_elements_file)
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, consensus_elements);
					}

					consensus_elements_file
						<< "#  Elements of consensus features extracted from " << in
						<< " on " << date_time_now << std::endl
						<< "# RT MZ Intensity Charge" << std::endl;
					for ( ConsensusMap::const_iterator cmit = consensus_map.begin(); cmit != consensus_map.end(); ++cmit )
					{
						consensus_elements_file << std::endl;
						for ( ConsensusFeature::const_iterator cfit = cmit->begin(); cfit != cmit->end(); ++cfit )
						{
							consensus_elements_file
								<< "H " << FeatureHandlePrinter(*cfit)
								<< "    " << ConsensusFeaturePrinter(*cmit) << std::endl;
						}
						// We repeat the first feature handle at the end of the list.
						// This way you can generate closed line drawings
						// See Gnuplot set datafile commentschars
						consensus_elements_file
							<< "L " << FeatureHandlePrinter(*cmit->begin())
							<< "    " << ConsensusFeaturePrinter(*cmit) << std::endl;
					}
					consensus_elements_file.close();
				}

				// ----------------------------------------------------------------------
			
				if ( !consensus_features.empty() )
				{
					std::ofstream consensus_features_file(consensus_features.c_str());
					if (!consensus_features_file)
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, consensus_features);
					}

					std::map<UInt,UInt> map_id_to_map_num;
					std::vector<UInt> map_num_to_map_id;
					std::vector<FeatureHandle> feature_handles;
					FeatureHandle feature_handle_NaN;
					feature_handle_NaN.setRT(std::numeric_limits<FeatureHandle::CoordinateType>::quiet_NaN());
					feature_handle_NaN.setMZ(std::numeric_limits<FeatureHandle::CoordinateType>::quiet_NaN());
					feature_handle_NaN.setIntensity(std::numeric_limits<FeatureHandle::IntensityType>::quiet_NaN());
					// feature_handle_NaN.setCharge(std::numeric_limits<Int>::max());

					for(ConsensusMap::FileDescriptions::const_iterator fdit = consensus_map.getFileDescriptions().begin();
							fdit != consensus_map.getFileDescriptions().end();
							++fdit
						 )
					{
						map_id_to_map_num[fdit->first] = map_num_to_map_id.size();
						map_num_to_map_id.push_back(fdit->first);
					}

					consensus_features_file
						<< "#  Consensus features extracted from " << in
						<< " on " << date_time_now << std::endl
						<< "# RT_cf MZ_cf Intensity_cf Charge_cf";
					for ( UInt fhindex = 0; fhindex < map_num_to_map_id.size(); ++fhindex )
					{
						const UInt map_id = map_num_to_map_id[fhindex];
						consensus_features_file
							<< "    RT_" << map_id
							<< " MZ_" << map_id
							<< " Intensity_" << map_id
							<< " Charge_" << map_id;
					}
					consensus_features_file << std::endl;

					for ( ConsensusMap::const_iterator cmit = consensus_map.begin(); cmit != consensus_map.end(); ++cmit )
					{
						{
							// please can anyone explain to me why putting the next two things into one statement doesnt work?
							std::vector<FeatureHandle> tmp(map_num_to_map_id.size(),feature_handle_NaN);
							feature_handles.swap(tmp);
						}
						consensus_features_file << ConsensusFeaturePrinter(*cmit);
						for ( ConsensusFeature::const_iterator cfit = cmit->begin(); cfit != cmit->end(); ++cfit )
						{
							feature_handles[map_id_to_map_num[cfit->getMapIndex()]] = *cfit;
						}
						for ( UInt fhindex = 0; fhindex < feature_handles.size(); ++fhindex )
						{
							consensus_features_file << "    " << FeatureHandlePrinter(feature_handles[fhindex]);
						}
						consensus_features_file << std::endl;
					}
					consensus_features_file.close();
				}

				return EXECUTION_OK;
      }
			else if (in_type == FileHandler::IDXML)
			{
				vector<ProteinIdentification> prot_ids;
				vector<PeptideIdentification> pep_ids;
				IdXMLFile().load(in, prot_ids, pep_ids);
				
				counter = 0;
				ofstream txt_out(out.c_str());

				for (vector<ProteinIdentification>::const_iterator it = prot_ids.begin(); it != prot_ids.end(); ++it)
				{
					String actual_id = it->getIdentifier();
					if (!getFlag_("peptides_only") && !getFlag_("peptides_only_csv"))
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
								if (!without_header_repetition)
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
								}
											
								// header of peptide hits
								if (without_header_repetition && counter == 0)
								{
									if (first_dim_rt)
									{
            				txt_out << "RT MZ Score Rank Sequence Charge AABefore AAAfter Accessions predicted_RT RT_first_dim predicted_RT_first_dim" << endl;
									}
									else
									{
            				txt_out << "RT MZ Score Rank Sequence Charge AABefore AAAfter Accessions predicted_RT" << endl;
            			}
            			++counter;
            		}
            		else if (counter == 0)
            		{
            			txt_out << "# Peptide Hits: Score, Rank, Sequence, Charge, AABefore, AAAfter, Accessions, predicted_RT" << endl;
            		}

            		for (vector<PeptideHit>::const_iterator ppit = pit->getHits().begin(); ppit != pit->getHits().end(); ++ppit)
            		{
            			if (without_header_repetition)
            			{
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
            			}
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
									if (ppit->metaValueExists("predicted_RT"))
									{
										txt_out << " " << ppit->getMetaValue("predicted_RT");
									}
									else
									{
										txt_out << " -1";
									}
									if (first_dim_rt)
									{
										if (pit->metaValueExists("first_dim_rt"))
										{
											txt_out << " " << pit->getMetaValue("first_dim_rt");
										}
										else
										{
											txt_out << " -1";
										}
										if (ppit->metaValueExists("predicted_RT_first_dim"))
										{
											txt_out << " " << ppit->getMetaValue("predicted_RT_first_dim");
										}
										else
										{
											txt_out << " -1";
										}										
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
}

int main( int argc, const char** argv )
{
	TOPPTextExporter t;
	return t.main(argc,argv);
}

/// @endcond

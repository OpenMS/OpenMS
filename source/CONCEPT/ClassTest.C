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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

namespace OpenMS
{
	namespace Internal
	{
		namespace ClassTest
		{
			bool all_tests = true;
			bool equal_files;
			bool newline = false;
			bool test = true;
			bool this_test;
			char line_buffer[65536];
			const char* version_string = 0;
			double absdiff = 0.;
			double absdiff_max = 0.;
			double absdiff_max_allowed = 1E-5;
			double ratio = 1.;
			double ratio_max = 1.;
			double ratio_max_allowed = 1. + 1E-5;
			int exception = 0;
			int line_num_1_max = -1;
			int line_num_2_max = -1;
			int start_section_line = 0;
			int test_count = 0;
			int test_line = 0;
			int verbose = 0;
			std::ifstream infile;
			std::ifstream templatefile;
			std::string add_message;
			std::string exception_message = "";
			std::string exception_name = "";
			std::string fuzzy_message;
			std::string test_name = "";
			std::vector<std::string> tmp_file_list;
			StringList whitelist;
		}
		
		namespace ClassTest
		{
			
			void setWhitelist(const char * const /* file */, const int line, const std::string& whitelist)
			{
				TEST::whitelist = StringList::create(whitelist);

				if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	
				{
					TEST::initialNewline();	
					std__cout << "    (line " << line <<															
						":  WHITELIST(\"" << whitelist << "\"):   whitelist is: " << TEST::whitelist << std::endl;																			
				}
				return;
			}

			void initialNewline()
			{
				if (!newline)
				{
					newline = true;
					std::cout << std::endl;
				}
				return;
			}
		
			void printWithPrefix(const std::string & text, const int marked)
			{
				std::istringstream is(text);
				std::string line;
				int line_number = 0;
				while ( std::getline(is,line) )
				{
					++line_number;
					std::cout << ( line_number == marked ? " # :|:  " : "   :|:  " ) << line << '\n';
				}
				return;
			}

			bool validate(const std::vector<std::string>& file_names)
			{
				bool passed = true;
				for (UInt i=0; i<file_names.size(); ++i)          								
				{																																									
					if (File::exists(file_names[i]))																
					{																																								
						switch(FileHandler::getType(file_names[i]))					
						{																																							
						case FileHandler::MZML:
							{																	
								if (!MzMLFile().isValid(file_names[i]))								
								{																																						
									std::cout << "Error: mzML file does not validate against XML schema '" << file_names[i] << "' - " << std::endl; 
									passed = false;																							
								}
								StringList errors, warnings;
								if (!MzMLFile().isSemanticallyValid(file_names[i], errors, warnings))								
								{																																						
									std::cout << "Error: mzML file semantically invalid '" << file_names[i] << "' - " << std::endl;
									for (UInt j=0; j<errors.size(); ++j)
									{
										std::cout << "Error - " << errors[j] << std::endl;
									}
									passed = false;																							
								}
							}
							break;	
						case FileHandler::MZDATA:																						
							if (!MzDataFile().isValid(file_names[i]))								
							{																																						
								std::cout << "Error: Invalid mzData file '" << file_names[i] << "' - " << std::endl; 
								passed = false;																							
							}
							break;																																			
						case FileHandler::MZXML:																											
							if (!MzXMLFile().isValid(file_names[i]))													
							{																																						
								std::cout << "Error: Invalid mzXML file '" << file_names[i] << "' - " << std::endl; 
								passed = false;																									
							}																																						
							break;																																			
						case FileHandler::FEATUREXML:																									
							if (!FeatureXMLFile().isValid(file_names[i]))													
							{																																						
								std::cout << "Error: Invalid FeatureXML file '" << file_names[i] << "' - " << std::endl; 
								passed = false;																									
							}																																						
							break;																																																																					
						case FileHandler::IDXML:																											
							if (!IdXMLFile().isValid(file_names[i]))													
							{																																						
								std::cout << "Error: Invalid IdXML file '" << file_names[i] << "' - " << std::endl; 
								passed = false;																									
							}																																						
							break;																																			
						case FileHandler::CONSENSUSXML:								
							if (!ConsensusXMLFile().isValid(file_names[i]))													
							{																																						
								std::cout << "Error: Invalid ConsensusXML file '" << file_names[i] << "' - " << std::endl; 
								passed = false;																									
							}																																						
							break;																																			
						case FileHandler::PARAM:																						
							if (!Param().isValid(file_names[i]))													
							{																																						
								std::cout << "Error: Invalid Param file '" << file_names[i] << "' - " << std::endl; 
								passed = false;																									
							}																																						
							break;																											
						case FileHandler::TRANSFORMATIONXML:																						
							if (!TransformationXMLFile().isValid(file_names[i]))													
							{																																						
								std::cout << "Error: Invalid TransformationXML file '" << file_names[i] << "' - " << std::endl; 
								passed = false;																									
							}																																						
							break;																																														
						default:																																			
							break;																																			
						}																																								
					}																																									
				}																												
				return passed;
			}
		
			std::string tmpFileName(const std::string& file, int line)
			{
				return String(file).prefix('.') + '_' + String(line) + ".tmp";
			}
		
 			bool isRealSimilar(double number_1, double number_2)
			{

				// Note: The original version of the stuff below was copied from
				// FuzzyStringComparator and then heavily modified for ClassTest.
				// Thus the case distinctions should be similar.

				absdiff = 0.;
				ratio = 0.;
				fuzzy_message.clear();

				// check if absolute difference is small
				absdiff = number_1 - number_2;
				if ( absdiff < 0 )
				{
					absdiff = -absdiff;
				}
				if ( absdiff > absdiff_max )
				{
					absdiff_max = absdiff;
				}
				// If absolute difference is small, large relative errors will be
				// tolerated in the cases below.  But a large absolute difference is
				// not an error, if relative error is small.  We do not jump out of
				// the case distinction here because we want to record the relative
				// error even in case of a successful comparison.
				bool is_absdiff_small = ( absdiff <= absdiff_max_allowed );
					
				if ( !number_1 )
				{ // number_1 is zero
					if (!number_2 )
					{ // both numbers are zero
						fuzzy_message = "both numbers are zero";
						return true;
					}
					else
					{
						if ( !is_absdiff_small )
						{
							fuzzy_message = "number_1 is zero, but number_2 is not small";
							return false;
						}
						else
						{
							fuzzy_message = "number_1 is zero, number_2 is small";
							return true;
						}
					}
				}
				else
				{ // number_1 is not zero
					if ( !number_2 )
					{
						if ( !is_absdiff_small )
						{
							fuzzy_message = "number_1 is not zero, but number_2 is";
							return false;
						}
						else
						{
							fuzzy_message = "number_2 is zero, but number_1 is not small";
							return true;
						}
					}
					else
					{ // both numbers are not zero
						ratio = number_1 / number_2;
						if ( ratio < 0. )
						{
							if ( !is_absdiff_small )
							{
								fuzzy_message = "numbers have different signs and difference is not small";
								return false;
							}
							else
							{
								fuzzy_message = "numbers have different signs, but difference is small";
								return true;
							}
						}
						else
						{ // ok, numbers have same sign, but we still need to check their ratio
							if ( ratio < 1. )
							{ // take reciprocal value
								ratio = 1. / ratio;
							}
							// by now, we are sure that ratio >= 1
							if ( ratio > ratio_max )
							{ // update running max
								ratio_max = ratio;
							}
							if ( ratio > ratio_max_allowed )
							{
								if ( !is_absdiff_small )
								{
									fuzzy_message = "ratio of numbers is large";
									return false;
								}
								else
								{
									fuzzy_message = "ratio of numbers is large, but numbers are small";
									return true;
								}
							}
							else
							{
								fuzzy_message = "ratio of numbers is small";
								return true;
							}
						}
					}
				}
				// We should never get here ... must have forgotten an condition branch above ... and then we need to fix that.
				fuzzy_message = "error: ClassTest.C:  You should never see this message.  Please report this bug along with the data that produced it.";
				return false;
			}

			bool isStringSimilar( const std::string & string_1, const std::string & string_2)
			{
				fuzzy_message.clear();
				FuzzyStringComparator fsc;
				fsc.setAcceptableAbsolute(absdiff_max_allowed);
				fsc.setAcceptableRelative(ratio_max_allowed);
				fsc.setVerboseLevel(2);
				fsc.setWhitelist(whitelist);
				std::ostringstream os;
				fsc.setLogDestination(os);
				fsc.use_prefix_ = true;

				bool result = fsc.compareStrings(string_1,string_2);

				fuzzy_message = os.str();
				absdiff = fsc.absdiff_max_;
				ratio = fsc.ratio_max_;
				line_num_1_max = fsc.line_num_1_max_;
				line_num_2_max = fsc.line_num_2_max_;

				return result;
			}

			bool isFileSimilar( const std::string & filename_1, const std::string & filename_2)
			{
				fuzzy_message.clear();
				FuzzyStringComparator fsc;
				fsc.setAcceptableAbsolute(absdiff_max_allowed);
				fsc.setAcceptableRelative(ratio_max_allowed);
				fsc.setVerboseLevel(2);
				fsc.setWhitelist(whitelist);
				std::ostringstream os;
				fsc.setLogDestination(os);
				fsc.use_prefix_ = true;

				bool result = fsc.compareFiles(filename_1,filename_2);

				fuzzy_message = os.str();
				absdiff = fsc.absdiff_max_;
				ratio = fsc.ratio_max_;
				line_num_1_max = fsc.line_num_1_max_;
				line_num_2_max = fsc.line_num_2_max_;

				return result;
			}

		}
	}
}

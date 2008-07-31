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

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

using namespace OpenMS;

namespace OpenMS
{
	namespace Internal
	{
		
		bool ClassTest::validate(const std::vector<std::string>& file_names)
		{
			bool passed = true;
			for (UInt i=0; i<file_names.size(); ++i)          								
			{																																									
				if (File::exists(file_names[i]))																
				{																																								
					switch(FileHandler::getType(file_names[i]))					
					{																																							
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
								std::cout << "Error: Invalid FeaturePairsXML file '" << file_names[i] << "' - " << std::endl; 
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
		
		
		std::string ClassTest::tmpFileName(const std::string& file, int line)
		{
			return String(file).prefix('.') + '_' + String(line) + ".tmp";
		}
		
	}
}

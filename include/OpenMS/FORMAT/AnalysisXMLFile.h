// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_ANALYSISXMLFILE_H
#define OPENMS_FORMAT_ANALYSISXMLFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ContactPerson.h>



#include <fstream>
#include <vector>

namespace OpenMS 
{
  /**
    @brief Used to load and store analysisXML files
    
    This class is used to load and store documents that implement 
    the schema of analysisXML files.
  
  	@ingroup FileIO
  */
  class AnalysisXMLFile
  {
    public:
      /// Constructor
      AnalysisXMLFile();
      /// Copy constructor
      AnalysisXMLFile(const AnalysisXMLFile& source);
      /// Destructor
      ~AnalysisXMLFile();
      
			/**
      	@brief Loads the identifications of an AnalysisXML file

				The information is read in and the information is stored in the
				corresponding variables
      */
      void load(const String& filename,
      					std::vector<ProteinIdentification>& protein_identifications, 
      					std::vector<Identification>& identifications, 
      					std::vector<float>& precursor_retention_times,
      					std::vector<float>& precursor_mz_values,
      					ContactPerson& contact_person)  	const throw (Exception::FileNotFound, 
  							 																							 Exception::ParseError);
      					 
			/**
      	@brief Loads the identifications of an AnalysisXML file

				The information is read in and the information is stored in the
				corresponding variables. This function offers the possibility to load 
				the predicted retention times and the predicted sigma that can be
				used to filter the peptides via the TOPP tool RTPredict.
      */
      void load(const String& filename, 
      					std::vector<ProteinIdentification>& protein_identifications, 
      					std::vector<Identification>& identifications, 
      					std::vector<float>& precursor_retention_times,
      					std::vector<float>& precursor_mz_values,
      					ContactPerson& contact_person,
      					std::map<String, double>& predicted_retention_times,
      					DoubleReal& predicted_sigma)  	const throw (Exception::FileNotFound, 
  							 																						 Exception::ParseError);
      					 
			/**
      	@brief Stores the data in an AnalysisXML file

				The data is read in and stored in the file 'filename'.
      */
      void store(String filename, 
	      				 const std::vector<ProteinIdentification>& protein_identifications, 
      					 const std::vector<Identification>& identifications, 
      					 const std::vector<float>& precursor_retention_times,
      					 const std::vector<float>& precursor_mz_values) const throw (Exception::UnableToCreateFile); 
			
			/**
      	@brief Stores the data in an AnalysisXML file

				The data is read in and stored in the file 'filename'.
      */
      void store(String filename, 
	      				 const std::vector<ProteinIdentification>& protein_identifications, 
      					 const std::vector<Identification>& identifications, 
      					 const std::vector<float>& precursor_retention_times,
      					 const std::vector<float>& precursor_mz_values,
      					 const ContactPerson& contact_person) const throw (Exception::UnableToCreateFile); 

			/**
      	@brief Stores the data in an AnalysisXML file

				The data is read in and stored in the file 'filename'.
      */
      void store(String filename, 
	      				 const std::vector<ProteinIdentification>& protein_identifications, 
      					 const std::vector<Identification>& identifications, 
      					 const std::vector<float>& precursor_retention_times,
      					 const std::vector<float>& precursor_mz_values,
      					 const ContactPerson& contact_person,
      					 const std::map<String, double>& predicted_retention_times,
      					 DoubleReal predicted_sigma) 
      	const throw (Exception::UnableToCreateFile); 

  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_ANALYSISXMLFILE_H

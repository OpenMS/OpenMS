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
// $Id: AnalysisXMLHandler.h,v 1.8 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_ANALYSISXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_ANALYSISXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>
#include <map>

namespace OpenMS
{
	namespace Internal
	{
  /**
    @brief Handler that is used for parsing AnalysisXML data
    
  */
  class AnalysisXMLHandler:
    public XMLHandler
  {
    public:
      /// Constructor
      AnalysisXMLHandler(std::vector<ProteinIdentification>& protein_identifications,
      									 std::vector<Identification>& identifications, 
      									 std::vector<float>& precursor_retention_times, 
      									 std::vector<float>& precursor_mz_values);
      /// Constructor
      AnalysisXMLHandler(const std::vector<ProteinIdentification>& protein_identifications,
      									 const std::vector<Identification>& identifications, 
      									 const std::vector<float>& precursor_retention_times, 
      									 const std::vector<float>& precursor_mz_values);
      /// Constructor
      AnalysisXMLHandler(std::vector<ProteinIdentification>* protein_identifications,
      									 std::vector<Identification>* identifications, 
      									 std::vector<float>* precursor_retention_times, 
      									 std::vector<float>* precursor_mz_values,
      									 ContactPerson* contact_person);
      /// Constructor
      AnalysisXMLHandler(std::vector<ProteinIdentification>* protein_identifications,
      									 std::vector<Identification>* identifications, 
      									 std::vector<float>* precursor_retention_times, 
      									 std::vector<float>* precursor_mz_values,
      									 ContactPerson* contact_person,
      									 std::map<String, double>* predicted_retention_times,
      									 DoubleReal* predicted_sigma);
      /// Constructor
      AnalysisXMLHandler(const std::vector<ProteinIdentification>& protein_identifications,
      									 const std::vector<Identification>& identifications, 
      									 const std::vector<float>& precursor_retention_times, 
      									 const std::vector<float>& precursor_mz_values,
      									 const ContactPerson& contact_person);
      /// Constructor
      AnalysisXMLHandler(const std::vector<ProteinIdentification>& protein_identifications,
      									 const std::vector<Identification>& identifications, 
      									 const std::vector<float>& precursor_retention_times, 
      									 const std::vector<float>& precursor_mz_values,
      									 const ContactPerson& contact_person,
      									 const std::map<String, double>& predicted_retention_times,
      									 DoubleReal predicted_sigma);
      /// Copy constructor
      AnalysisXMLHandler(const AnalysisXMLHandler& source);
      /// Destructor
      ~AnalysisXMLHandler();
      
      /// Assignment operator
      AnalysisXMLHandler& operator = (const AnalysisXMLHandler& source);
      
      /// Equality operator
      bool operator == (const AnalysisXMLHandler& source) const;
      
      /// Equality operator
      bool operator != (const AnalysisXMLHandler& source) const;
      
      void writeTo(std::ostream& os);
      	
      virtual bool startElement(const QString & uri, const QString & local_name,
												const QString & qname, const QXmlAttributes & attributes );
			///
      virtual bool endElement( const QString & uri, const QString & local_name,
											 const QString & qname ); 
						
		  bool characters( const QString & chars );
    protected:
      std::vector<ProteinIdentification>* protein_identifications_;
      std::vector<Identification>* identifications_;
      std::vector<float>* precursor_retention_times_;
      std::vector<float>* precursor_mz_values_;
      ProteinHit actual_protein_hit_;
      std::vector<ProteinHit> actual_protein_hits_;
      PeptideHit actual_peptide_hit_;
      std::vector<PeptideHit> actual_peptide_hits_;
			UnsignedInt peptide_identification_index_;
			UnsignedInt protein_identification_index_;
			bool inside_peptide_;
			ContactPerson* contact_person_;      	
      const std::vector<ProteinIdentification> const_protein_identifications_;
      const std::vector<Identification> const_identifications_;
      const std::vector<float> const_precursor_retention_times_;
      const std::vector<float> const_precursor_mz_values_;
			const ContactPerson const_contact_person_;
			const std::map<String, double> const_predicted_retention_times_;
			String tag_;      	
			UnsignedInt charge_identification_index_;
			bool inside_protein_;
			bool inside_global_protein_;
			std::vector<UnsignedInt> actual_peptide_indices_;
			std::map<String, double>* predicted_retention_times_;
			DoubleReal* predicted_sigma_;
			DoubleReal const_predicted_sigma_;
			std::vector< String > date_times_temp_;
			UnsignedInt date_times_counter_;
			
		private:
		    UnsignedInt getDateGroupIndex(DateTime 												date_time,
			  															std::map< String , UnsignedInt> date_times);
				void writePeptideHit(std::ostream& os, 
		 													String shift,
		 													PeptideHit hit,
		 													Real significance_threshold,
		 													UnsignedInt identification_index,
		 													SignedInt charge, 
		 													Real precursor_retention_time,
		 													Real precursor_mz,
		 													DateTime date_time,
					  									std::map< String , UnsignedInt> date_times);


  };

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_ANALYSISXMLHANDLER_H

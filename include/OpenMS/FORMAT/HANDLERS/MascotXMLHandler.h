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
// $Id: MascotXMLHandler.h,v 1.1 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MASCOTXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MASCOTXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/METADATA/Identification.h>

#include <vector>

namespace OpenMS
{
	namespace Internal
	{
  /**
    @brief Handler that is used for parsing MascotXML data
    
  */
  class MascotXMLHandler:
    public XMLHandler
  {
    public:
      /// Constructor
      MascotXMLHandler(ProteinIdentification*       protein_identification,
      								 std::vector<Identification>* identifications, 
      								 std::vector<float>* 					precursor_retention_times, 
      								 std::vector<float>* 					precursor_mz_values);

      /// Copy constructor
      MascotXMLHandler(const MascotXMLHandler& source);

      /// Destructor
      ~MascotXMLHandler();
      
      /// Assignment operator
      MascotXMLHandler& operator = (const MascotXMLHandler& source);
      
      /// Equality operator
      bool operator == (const MascotXMLHandler& source) const;
      
      /// Equality operator
      bool operator != (const MascotXMLHandler& source) const;
      
      virtual bool startElement(const QString & uri, const QString & local_name,
												const QString & qname, const QXmlAttributes & attributes );
			///
      virtual bool endElement( const QString & uri, const QString & local_name,
											 const QString & qname ); 
						
		  bool characters( const QString & chars );
    protected:
    	ProteinIdentification* protein_identification_;
      std::vector<Identification>* identifications_;
      std::vector<float>* precursor_retention_times_;
      std::vector<float>* precursor_mz_values_;
      ProteinHit actual_protein_hit_;
      std::vector<ProteinHit> actual_protein_hits_;
      PeptideHit actual_peptide_hit_;
      std::vector<PeptideHit> actual_peptide_hits_;
			UnsignedInt peptide_identification_index_;
			UnsignedInt protein_identification_index_;
    	const ProteinIdentification const_protein_identification_;
      const std::vector<Identification> const_identifications_;
      const std::vector<float> const_precursor_retention_times_;
      const std::vector<float> const_precursor_mz_values_;
			String tag_;      	
			bool inside_protein_;
			std::vector<UnsignedInt> actual_peptide_indices_;
			DateTime date_;
			String date_time_string_;
  };

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_MASCOTXMLHANDLER_H

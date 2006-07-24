/// -*- Mode: C++; tab-width: 2; -*-
/// vi: set ts=2:
///
/// --------------------------------------------------------------------------
///                   OpenMS Mass Spectrometry Framework
/// --------------------------------------------------------------------------
///  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
///
///  This library is free software; you can redistribute it and/or
///  modify it under the terms of the GNU Lesser General Public
///  License as published by the Free Software Foundation; either
///  version 2.1 of the License, or (at your option) any later version.
///
///  This library is distributed in the hope that it will be useful,
///  but WITHOUT ANY WARRANTY; without even the implied warranty of
///  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
///  Lesser General Public License for more details.
///
///  You should have received a copy of the GNU Lesser General Public
///  License along with this library; if not, write to the Free Software
///  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
///
/// --------------------------------------------------------------------------
/// $Id: Outfile.h,v 1.0 2006/07/12 15:58:59 martinlangwisch Exp $
/// $Author: martinlangwisch $
/// $Maintainer: Martin Langwisch $
/// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_OUTFILE_H
#define OPENMS_FORMAT_OUTFILE_H


#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <cstring>

namespace OpenMS
{
  /**
    @brief Representation of an  outfile
    
    This class serves to read in an  outfile. The information can be 
    retrieved via the >> operator. 
  
  	@ingroup FileIO
  */
  class Outfile
  {
    public:
			/// empty constructor
      Outfile();
			
      /// Copy constructor
      Outfile(const Outfile& Outfile);
			
			/// Destructor
      ~Outfile();
			
      /// true if the search was successfull, false otherwise
      bool ok() const;
			
      /// fills a Identification object
      Outfile& operator>>(Identification& Identification);
			
      /// fills a PeptideHit object
      Outfile& operator>>(PeptideHit& peptide_hit);
			
      /// fills a ProteinHit object
      Outfile& operator>>(ProteinHit& protein_hit);
			
			/// Assignment operator
	    Outfile& operator=(const Outfile& source);
			
      /// returns the retention time of the  search
      const std::vector<float>& getPrecursorRetentionTimes() const;
			
      /// sets the retention time of the  search
      void setPrecursorRetentionTimes(const std::vector<float>& precursor_retention_times);      
			
      /// returns the m/z of the precursor peak of the  search
      const std::vector<float>& getPrecursorMZValues() const;
			
      /// sets the m/z of the precursor peak of the  search
      void setPrecursorMZValues(const std::vector<float>& mz);
			
      /// returns the Identification instances of the  search
      const std::vector<Identification>& getIdentifications() const;
			
      /// sets the Identification instances of the  search
      void setIdentifications(const std::vector<Identification>& queries);
			
      /// returns the ProteinIdentification instances of the  search
      const ProteinIdentification& getProteinIdentification() const;
			
      /// sets the ProteinIdentification instances of the  search
      void setProteinIdentification(const ProteinIdentification& protein_ids);
			
    protected:			
			/// get the accession and accession type
			void get_ac_and_ac_type(String line, const std::string& filename, std::string& accession, std::string& accession_type) throw (Exception::ParseError);
			
			/// given a vector of peptide hits, either insert the new peptide hit or update its ProteinHits, returns whether an update took place
			bool updatePeptideHits(PeptideHit& peptide_hit, std::vector< PeptideHit >& peptide_hits);
			
			/// the Identification information
			std::vector< Identification > queries_;
			
			/// the ProteinIdentification information
			ProteinIdentification protein_ids_;
			
			/// list of the peptide hits (sorted by score)
			std::vector< PeptideHit > peptide_hits_;
			
			/// list of the protein hits (sorted by score)
			std::vector< ProteinHit > protein_hits_;
			
			/// the retention time
			std::vector< float > precursor_retention_times_;
			
      /// the mass of the precursor
      std::vector< float > precursor_mz_values_;
      
      /// flag that states if the search worked
      bool ok_;
			
			/// iterator pointing to the current db search
			std::vector< Identification >::iterator curr_query_;
			
			/// iterator pointing to the current hit
			std::vector< PeptideHit >::iterator curr_peptide_hit_;
			
			/// iterator pointing to the current hit
			std::vector< ProteinHit >::iterator curr_protein_hit_;
			
   };
	
} //namespace OpenMS

#endif // OPENMS_FORMAT_OUTFILE_H

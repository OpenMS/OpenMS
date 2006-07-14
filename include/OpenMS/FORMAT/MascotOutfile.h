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
// $Maintainer: Nico Pfeiffer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MASCOTOUTFILE_H
#define OPENMS_FORMAT_MASCOTOUTFILE_H

#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <fstream>
#include <iostream>
#include <string>
#include <map>

namespace OpenMS
{
  /**
    @brief Representation of a Mascot outfile
    
    This class serves to read in a Mascot outfile. The information can be 
    retrieved via the >> operator. 
  	
  	@todo adapt to common interface: load/store (Nico)
  	
  	@ingroup FileIO
  */
  class MascotOutfile
  {
    public:
      
      /// Constructor
      MascotOutfile(const std::string& filename, Real p = 0.05) 
      	throw (Exception::ParseError);

      /// Copy constructor
      MascotOutfile(const MascotOutfile& mascotoutfile);

      /// Destructor
      ~MascotOutfile();

      /// true if the search was successfull, false otherwise
      bool ok() const;

      /// fills a Identification object
      MascotOutfile& operator>>(Identification& identification);

      /// fills a PeptideHit object
      MascotOutfile& operator>>(PeptideHit& peptide_hit);

      /// fills a ProteinHit object
      MascotOutfile& operator>>(ProteinHit& protein_hit);

			/// Assignment operator
	    MascotOutfile& operator=(const MascotOutfile& source);
		
      /// returns the retention time of the Mascot search
      const std::vector<float>& getPrecursorRetentionTimes() const;

      /// sets the retention time of the Mascot search
      void setPrecursorRetentionTimes(const std::vector<float>& precursor_retention_times);      

      /// returns the m/z of the precursor peak of the Mascot search
      const std::vector<float>& getPrecursorMZValues() const;

      /// sets the m/z of the precursor peak of the Mascot search
      void setPrecursorMZValues(const std::vector<float>& mz);      

      /// returns the Identification instances of the Mascot search
      const std::vector<Identification>& getIdentifications() const;

      /// sets the Identification instances of the Mascot search
      void setIdentifications(const std::vector<Identification>& identifications);      

    protected:

			/// the identification information
			std::vector<Identification> db_searches_;

			/// list of the peptide hits (sorted by score)
			std::vector<PeptideHit> peptide_hits_;

			/// list of the protein hits (sorted by score)
			std::vector<ProteinHit> protein_hits_;

			/// the retention time
			std::vector<float> precursor_retention_times_;

			/// iterator pointing to the current hit
			std::vector<PeptideHit>::iterator curr_peptide_hit_;
				
			/// iterator pointing to the current hit
			std::vector<ProteinHit>::iterator curr_protein_hit_;
				
      /// the mass of the precursor
      std::vector<float> precursor_mz_values_;
      
      /// flag that states if the search worked
      bool ok_;
      
   };

} //namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTOUTFILE_H

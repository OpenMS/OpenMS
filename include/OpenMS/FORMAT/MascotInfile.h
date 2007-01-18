// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FORMAT_MASCOTINFILE_H
#define OPENMS_FORMAT_MASCOTINFILE_H

#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <string>
#include <vector>

namespace OpenMS
{
	/**
		@brief Mascot input file adapter.
		
		Creates a file that can be used for Mascot search from a peak list or a whole experiment.
  	  	
  	@ingroup FileIO
	*/
  class MascotInfile
  {
    public:

			/// constructor
			MascotInfile();

			/// constructor
			~MascotInfile();

			/// stores the peak list in a MascotInfile that can be used as input for MASCOT shell execution
			void store(const std::string& filename, const DPeakArray<1>& spec, double mz , double retention_time, std::string search_title);		

			/// stores the experiment data in a MascotInfile that can be used as input for MASCOT shell execution
			void store(const std::string& filename,
								 const MSExperiment< DPeak<1> >& experiment, 
								 std::string search_title);
														
			/// returns the boundary used for the MIME format
			const std::string& getBoundary();
		  /// sets the boundary used for the MIME format.<br>By default a 22 character random string is used
		  void setBoundary(const std::string& boundary);

			/// returns the DB to use
			const std::string& getDB();
		  /// sets the DB to use (default: MSDB). See &lt;mascot path&gt;/config/mascot.dat in "Databases" section for possible settings
		  void setDB(const std::string& db);

			/// returns the search type
		  const std::string& getSearchType();
		  /// sets the seach type (default: MIS). So far only MIS is supported!<br>Valid types are "MIS" (MS/MS Ion Search), "PMF" (Peptide Mass Fingerprint) , "SQ" (Sequence Query)
		  void setSearchType(const std::string& search_type);

		  /// returns the number of hits to report back
		  const std::string& getHits();
		  /// sets the number of hits to report back (default: 20)
		  void setHits(const std::string& hits);

			/// returns the enzyme used for cleavage
		  const std::string& getCleavage();
		  /// sets the enzyme used for cleavage (default: Trypsin). <BR>See &lt;mascot path&gt;/config/enzymes for possible settings.
		  void setCleavage(const std::string& cleavage);

			/// returns the used mass type ("Monoisotopic" or "Average")
		  const std::string& getMassType();
		  /// sets the used mass type "Monoisotopic" or "Average" (default: Monoisotopic)
		  void setMassType(const std::string& mass_type);

			/// returns a vector containing the fixed modifications (default: none)
		  const std::vector<String>& getModifications();
		  /// sets the fixed modifications (default: none). <BR>See &lt;mascot path&gt;/config/mod_file for possible settings.
		  void setModifications(const std::vector<String>& mods);

			/// returns a vector containing the variable modifications (default: none)
		  const std::vector<String>& getVariableModifications();
		  /// sets the fixed modifications (default: none). <BR>See &lt;mascot path&gt;/config/mod_file for possible settings.
		  void setVariableModifications(const std::vector<String>& mods);

			/// returns the instrument type
		  const std::string& getInstrument();
		  /// sets the instrument type (Default: Default). <BR>Possible instruments: ESI-QUAD-TOF, MALDI-TOF-PSD, ESI-TRAP, ESI-QUAD, ESI-FTICR, MALDI-TOF-TOF, ESI-4SECTOR, FTMS-ECD, MALDI-QUAD-TOF, MALDI-QIT-TOF
		  void setInstrument(const std::string& instrument);

			/// returns the number of allowed missed cleavages
		  UnsignedInt getMissedCleavages();
		  /// sets the number of allowed missed cleavages (default: 1)
		  void setMissedCleavages(UnsignedInt missed_cleavages);

			/// returns the precursor mass tolerance
		  float getPrecursorMassTolerance();
		  /// sets the precursor mass tolerance in Da (default: 2.0)
		  void setPrecursorMassTolerance(float precursor_mass_tolerance);

			/// returns the peak mass tolerance in Da
		  float getPeakMassTolerance();
		  /// sets the peak mass tolerance in Da (default: 1.0)
		  void setPeakMassTolerance(float ion_mass_tolerance);

			/// returns the taxonomy
		  const std::string& getTaxonomy();
		  /// sets the taxonomy (default: All entries). <BR>See &lt;mascot path&gt;/config/taxonomy for possible settings.
		  void setTaxonomy(const std::string& taxonomy);

		  /// returns the Mascot form version
			const std::string& getFormVersion();
		  /// sets the Mascot form version (default: 1.01)
		  void setFormVersion(const std::string& form_version);

		  /// returns the charges
			const std::string& getCharges();
		  /// sets the charges (default: 1+, 2+ and 3+)
		  void setCharges(std::vector<SignedInt>& charges);

    protected:
			/// parent mass
			double mz_;

			/// charge states to use
			String charges_;

			/// the search title of the mascot search
			std::string search_title_;

			/// the DB to search in
			std::string db_;

			/// search type: MIS, SQ or PMF
			std::string search_type_;

			/// number of hits to report
			std::string hits_;

			/// Enzyme used for cleavage
			std::string cleavage_;

			/// Monoisotopic/average mass
			std::string mass_type_;

			/// fixed Modifications
			std::vector<String> mods_;

			/// variable Modifications
			std::vector<String> variable_mods_;

			/// the used instument
			std::string instrument_;

			/// number of missed cleavages
			UnsignedInt missed_cleavages_;

			/// precursor mass toerance in Da
			float precursor_mass_tolerance_;

			/// m/z tolerance of ions  in Da
			float ion_mass_tolerance_;

			/// taxonomy
			std::string taxonomy_;

			/// form version
			std::string form_version_;

			/// the boundary used for the MIME format
			std::string boundary_;

			/// the retention time
			double retention_time_;

			/// writes a parameter header
			void writeParameterHeader_(const std::string& name, FILE* fp, bool line_break = true);

			/// writes the full header
			void writeHeader_(FILE* fp);
			
			/// writes the spectrum
			void writeSpectrum_(FILE* fp,
													const std::string& filename,
													const DPeakArray<1>& peaks);
						
			/// writes the MSExperiment
			void writeMSExperiment_(FILE* fp, 
															const std::string& filename, 
															const MSExperiment< DPeak<1> >& experiment);

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTINFILE_H

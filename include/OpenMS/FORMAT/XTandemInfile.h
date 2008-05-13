// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XTANDEMINFILE_H
#define OPENMS_FORMAT_XTANDEMINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/HANDLERS/XTandemInfileXMLHandler.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>

#include <vector>
#include <fstream>

namespace OpenMS
{
	/**
		@brief XTandem input file adapter
		
		This class is able to create a X!Tandem configuration file for a search
  	  	
  	@ingroup FileIO
	*/
  class XTandemInfile
  {
    public:

			enum ERROR_UNIT
			{
				DALTONS = 0,
				PPM
			};
						
			enum MASS_TYPE
			{
				MONOISOTOPIC = 0,
				AVERAGE
			};
		

			/// constructor
			XTandemInfile();

			/// constructor
			virtual ~XTandemInfile();

			//<note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>




			//<note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
			//<note type="input" label="spectrum, parent monoisotopic mass error plus">100</note>
			//<note type="input" label="spectrum, parent monoisotopic mass error minus">100</note>
			//<note type="input" label="spectrum, parent monoisotopic mass isotope error">yes</note>
			//<note type="input" label="spectrum, fragment monoisotopic mass error units">Daltons</note>
			//<note>The value for this parameter may be 'Daltons' or 'ppm': all other values are ignored</note>
			//<note type="input" label="spectrum, parent monoisotopic mass error units">ppm</note>
			//<note>The value for this parameter may be 'Daltons' or 'ppm': all other values are ignored</note>
			//<note type="input" label="spectrum, fragment mass type">monoisotopic</note>
			//<note>values are monoisotopic|average </note>
			void setFragmentMassTolerance(double tolerance);

			double getFragmentMassTolerance() const;

			void setPrecursorMassTolerancePlus(double tol);

			double getPrecursorMassTolerancePlus() const;

			void setPrecursorMassToleranceMinus(double tol);

			double getPrecursorMassToleranceMinus() const;

			void setPrecursorMonoisotopicError(MASS_TYPE mono_isotopic);

			MASS_TYPE getPrecursorMonoisotopicError() const;

			void setFragmentMassErrorUnit(ERROR_UNIT unit);

			ERROR_UNIT getFragmentMassErrorUnit() const;

			void setPrecursorMassErrorUnit(ERROR_UNIT unit);

			ERROR_UNIT getPrecursorMassErrorUnit() const;
	
			void setNumberOfThreads(UInt threads);

			UInt getNumberOfThreads() const;
			
			void setFixedModifications(const ModificationDefinitionsSet& mods);

			const ModificationDefinitionsSet& getFixedModifications() const;

			void setVariableModifications(const ModificationDefinitionsSet& mods);

			const ModificationDefinitionsSet& getVariableModifications() const;
			
			void setOutputFilename(const String& output);

			const String& getOutputFilename() const;

			void setInputFilename(const String& input_file);

			const String& getInputFilename() const;

			void setTaxonomyFilename(const String& filename);

			const String& getTaxonomyFilename() const;

			void setDefaultParametersFilename(const String& filename);

			const String& getDefaultParametersFilename() const;

			void setTaxon(const String& taxon);

			const String& getTaxon() const;

			void setMaxPrecursorCharge(Int max_charge);

			Int getMaxPrecursorCharge() const;
			
			void setNumberOfMissedCleavages(UInt missed_cleavages);

			UInt getNumberOfMissedCleavages() const;

			void setMaxValidEValue(double value);

			double getMaxValidEValue() const;
			
			void write(const String& filename) throw (Exception::UnableToCreateFile);

			void load(const String& filename) throw (Exception::FileNotFound, Exception::ParseError);

    protected:


			void writeTo_(std::ostream& os);

			void writeNote_(std::ostream& os, const String& type, const String& label, const String& value);

			void writeNote_(std::ostream& os, const String& type, const String& label, const char* value);

			void writeNote_(std::ostream& os, const String& type, const String& label, bool value);

			double fragment_mass_tolerance_;

			double precursor_mass_tolerance_plus_;

			double precursor_mass_tolerance_minus_;

			MASS_TYPE precursor_monoisotopic_error_;

			ERROR_UNIT precursor_mass_error_unit_;

			ERROR_UNIT fragment_mass_error_unit_;

			MASS_TYPE fragment_mass_type_;

			UInt max_precursor_charge_;
			
			double precursor_lower_mz_;

			double fragment_lower_mz_;

			UInt number_of_threads_;

			UInt batch_size_;
			
			ModificationDefinitionsSet fixed_modifications_;

			ModificationDefinitionsSet variable_modifications_;

			String input_filename_;

			String output_filename_;

			String taxonomy_file_;
		
			String taxon_;

			String cleavage_site_;

			// refinement
			bool refine_;

			double refine_max_valid_evalue_;


			// scoring
			UInt number_of_missed_cleavages_;

			String default_parameters_file_;
			
			// output parameters
			double max_valid_evalue_;
		

			//<note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
			std::vector<Internal::XTandemInfileNote> notes_;	
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTINFILE_H

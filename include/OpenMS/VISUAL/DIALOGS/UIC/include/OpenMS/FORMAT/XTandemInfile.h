// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XTANDEMINFILE_H
#define OPENMS_FORMAT_XTANDEMINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/HANDLERS/XTandemInfileXMLHandler.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/FORMAT/XMLFile.h>

namespace OpenMS
{
	/**
		@brief XTandem input file adapter
		
		This class is able to create a X!Tandem configuration file for a search
  	  	
  	@ingroup FileIO
	*/
  class OPENMS_DLLAPI XTandemInfile
  	: public Internal::XMLFile
  {
    public:

			/// error unit, either Da or ppm
			enum ErrorUnit
			{
				DALTONS = 0,
				PPM
			};
			
			/// Mass type of the precursor, either monoisotopic or average
			enum MassType
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

			/// setter for the fragment mass tolerance
			void setFragmentMassTolerance(double tolerance);

			/// returns the fragment mass tolerance
			double getFragmentMassTolerance() const;

			/// sets the precursor mass tolerance (plus only)
			void setPrecursorMassTolerancePlus(double tol);

			/// returns the precursor mass tolerance (plus only)
			double getPrecursorMassTolerancePlus() const;

			/// set the precursor mass tolerance (minus only)
			void setPrecursorMassToleranceMinus(double tol);

			/// returns the precursor mass tolerance (minus only)
			double getPrecursorMassToleranceMinus() const;

			/// sets the precursor mass type
			void setPrecursorErrorType(MassType mono_isotopic);

			/// returns the precursor mass type
			MassType getPrecursorErrorType() const;

			/// sets the fragment mass error unit (Da, ppm)
			void setFragmentMassErrorUnit(ErrorUnit unit);

			/// returns the fragment mass error unit (Da, ppm)
			ErrorUnit getFragmentMassErrorUnit() const;

			/// sets the precursor mass error unit (Da, ppm)
			void setPrecursorMassErrorUnit(ErrorUnit unit);

			/// returns the precursor mass error unit (Da, ppm)
			ErrorUnit getPrecursorMassErrorUnit() const;

			/// sets the number of threads used during the identifications
			void setNumberOfThreads(UInt threads);

			/// returns the number of threads
			UInt getNumberOfThreads() const;
			
			/// sets the modifications using a modification definitions set
			void setModifications(const ModificationDefinitionsSet& mods);

			/// returns the modifications set, using a modification definitions set
			const ModificationDefinitionsSet& getModifications() const;

			/// sets the output filename
			void setOutputFilename(const String& output);

			/// returns the output filename
			const String& getOutputFilename() const;

			/// sets the input filename
			void setInputFilename(const String& input_file);

			/// returns the input filename
			const String& getInputFilename() const;

			/// set the filename of the taxonomy file
			void setTaxonomyFilename(const String& filename);

			/// returns the filename of the taxonomy file
			const String& getTaxonomyFilename() const;

			/// sets the default paramters file
			void setDefaultParametersFilename(const String& filename);

			/// returns the default parameters file
			const String& getDefaultParametersFilename() const;

			/// sets the taxon used in the taxonomy file
			void setTaxon(const String& taxon);

			/// returns the taxon used in the taxonomy file
			const String& getTaxon() const;

			/// sets the max precursor charge 
			void setMaxPrecursorCharge(Int max_charge);

			/// returns the max precursor charge
			Int getMaxPrecursorCharge() const;
			
			/// sets the number of missed cleavages allowed
			void setNumberOfMissedCleavages(UInt missed_cleavages);

			/// returns the number of missed cleavages allowed
			UInt getNumberOfMissedCleavages() const;

			/// sets the max valid E-value allowed in the list
			void setMaxValidEValue(double value);

			/// returns the max valid E-value allowed in the list
			double getMaxValidEValue() const;
			
			/** writes the XTandemInfile to the given file

					@param filename the name of the file which is written
					@throw UnableToCreateFile is thrown if the given file could not be created
			*/
			void write(const String& filename);

			/** read the information from the given filename
	
					@param filename the file which should be read from
					@throw FileNotFound is thrown if the given file could not be found
					@throw ParseError is thrown if the given file could not be parsed
			*/
			void load(const String& filename);

    protected:

			XTandemInfile(const XTandemInfile& rhs);

			XTandemInfile& operator = (const XTandemInfile& rhs);

			void writeTo_(std::ostream& os);

			void writeNote_(std::ostream& os, const String& type, const String& label, const String& value);

			void writeNote_(std::ostream& os, const String& type, const String& label, const char* value);

			void writeNote_(std::ostream& os, const String& type, const String& label, bool value);

			double fragment_mass_tolerance_;

			double precursor_mass_tolerance_plus_;

			double precursor_mass_tolerance_minus_;

			MassType precursor_mass_type_;

			ErrorUnit precursor_mass_error_unit_;

			ErrorUnit fragment_mass_error_unit_;

			MassType fragment_mass_type_;

			UInt max_precursor_charge_;
			
			double precursor_lower_mz_;

			double fragment_lower_mz_;

			UInt number_of_threads_;

			UInt batch_size_;
			
			ModificationDefinitionsSet modifications_;

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

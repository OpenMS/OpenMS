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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SPECTRUMSETTINGS_H
#define OPENMS_METADATA_SPECTRUMSETTINGS_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Identification.h>

#include <map>
#include <vector>

namespace OpenMS 
{
	/**
		@brief Representation of 1D spectrum settings. 
		
		It contains the metadata about spectrum specific instrument settings,
		acquisition settings, description of the meta values used in the peaks and precursor info.
		
		Precursor info should only be used if this spectrum is a tandem-MS spectrum.
		The precursor spectrum is the first spectrum before this spectrum, that has a lower MS-level than
		the current spectrum.
		
		@ingroup Metadata
	*/
  class SpectrumSettings 
  {
    public:
    	/**
    		@brief spectrum type
    		
    		If set to unknown see SpectrumType of the ProcessingMethod ( in MSExperiment ).
    	*/
    	enum SpectrumType {UNKNOWN, PEAKS, RAWDATA,SIZE_OF_SPECTRUMTYPE};
			/// Names of spectrum types
			static const std::string NamesOfSpectrumType[SIZE_OF_SPECTRUMTYPE];
   	
    	/// Constructor
      SpectrumSettings();
      /// Copy constructor
      SpectrumSettings(const SpectrumSettings& source);
      /// Destructor
      ~SpectrumSettings();
      
      // Assignment operator
      SpectrumSettings& operator= (const SpectrumSettings& source);

      /// Equality operator
      bool operator== (const SpectrumSettings& rhs) const;
      /// Equality operator
      bool operator!= (const SpectrumSettings& rhs) const;

			/**
				@brief returns the spectrum type
				
				If the type is 'UNKNOWN', a general type for all spectra of an experiment might be stored in the ProcessingMethod instance of an ExperimentalSettings .
			*/
      SpectrumType getType() const;
      /**
				@brief sets the spectrum type
				
				If the type is 'UNKNOWN', a general type for all spectra of an experiment might be stored in the ProcessingMethod instance of an ExperimentalSettings .
			*/
      void setType(SpectrumType type);

			/// returns the free-text comment
      const String& getComment() const;
      /// sets the free-text comment
      void setComment(const String& comment);
			
			/// returns a const reference to the instrument settings of the current spectrum
      const InstrumentSettings& getInstrumentSettings() const;
      /// returns a mutable reference to the instrument settings of the current spectrum
      InstrumentSettings& getInstrumentSettings();
      /// sets the instrument settings of the current spectrum
      void setInstrumentSettings(const InstrumentSettings& instrument_settings);
			
			/// returns a const reference to the acquisition info
      const AcquisitionInfo& getAcquisitionInfo() const;
      /// returns a mutable reference to the acquisition info
      AcquisitionInfo& getAcquisitionInfo();
      /// sets the acquisition info
      void setAcquisitionInfo(const AcquisitionInfo& acquisition_info);
			
			/// returns a const reference to the source file
      const SourceFile& getSourceFile() const;
      /// returns a mutable reference to the source file
      SourceFile& getSourceFile();
      /// sets the source file
      void setSourceFile(const SourceFile& source_file);
			
			/// returns a const reference to the description of the meta values used in the peaks
      const std::map<String,MetaInfoDescription>& getMetaInfoDescriptions() const;
      /// returns a mutable reference to the description of the meta values used in the peaks
      std::map<String,MetaInfoDescription>& getMetaInfoDescriptions();
      /// sets the description of the meta values used in the peaks
      void setMetaInfoDescriptions(const std::map<String,MetaInfoDescription>& meta_info_descriptions);

			/// returns a const reference to the precursor
      const Precursor& getPrecursor() const;
      /// returns a mutable reference to the precursor
      Precursor& getPrecursor();
      /// sets the precursor
      void setPrecursor(const Precursor& precursor);
			
      /// returns a const reference to the Identification vector
	    const std::vector<Identification>& getIdentifications() const;	    	
      /// returns a mutable reference to the Identification vector
	    std::vector<Identification>& getIdentifications();
	    /// sets the Identification vector
	    void setIdentifications(const std::vector<Identification>& identifications);	

    protected:
    	SpectrumType type_;
      String comment_;
      InstrumentSettings instrument_settings_;
      SourceFile source_file_;
      AcquisitionInfo acquisition_info_;
      std::map<String,MetaInfoDescription> meta_info_descriptions_;
      Precursor precursor_;
	    std::vector<Identification> identification_;
  };

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const SpectrumSettings& spec);

} // namespace OpenMS

#endif // OPENMS_METADATA_SPECTRUMSETTINGS_H

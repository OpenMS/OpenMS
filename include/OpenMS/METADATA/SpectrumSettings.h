// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SPECTRUMSETTINGS_H
#define OPENMS_METADATA_SPECTRUMSETTINGS_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Product.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/DataProcessing.h>

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
  class OPENMS_DLLAPI SpectrumSettings
  	: public MetaInfoInterface
  {
  	
    public:
    
    	///Spectrum peak type
    	enum SpectrumType
    	{
    		UNKNOWN,	///< Unknown spectrum type
    		PEAKS,		///< Peak data (also called centroided data or stick data)
    		RAWDATA,	///< Raw data (also called profile data)
    		SIZE_OF_SPECTRUMTYPE
    	};
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

      /// merge another spectrum setting into this one (data is usually appended, except for spectrum type which needs to be unambiguous to be kept)
      void unify(const SpectrumSettings& rhs);

			///returns the spectrum type
      SpectrumType getType() const;
      ///sets the spectrum type
      void setType(SpectrumType type);

			/// returns the native identifier for the spectrum, used by the acquisition software.
      const String& getNativeID() const;
      /// sets the native identifier for the spectrum, used by the acquisition software.
      void setNativeID(const String& native_id);

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
			
			/// returns a const reference to the precursors
      const std::vector<Precursor>& getPrecursors() const;
      /// returns a mutable reference to the precursors
      std::vector<Precursor>& getPrecursors();
      /// sets the precursors
      void setPrecursors(const std::vector<Precursor>& precursors);

			/// returns a const reference to the products
      const std::vector<Product>& getProducts() const;
      /// returns a mutable reference to the products
      std::vector<Product>& getProducts();
      /// sets the products
      void setProducts(const std::vector<Product>& products);
			
      /// returns a const reference to the PeptideIdentification vector
	    const std::vector<PeptideIdentification>& getPeptideIdentifications() const;	    	
      /// returns a mutable reference to the PeptideIdentification vector
	    std::vector<PeptideIdentification>& getPeptideIdentifications();
	    /// sets the PeptideIdentification vector
	    void setPeptideIdentifications(const std::vector<PeptideIdentification>& identifications);	

			/// returns a const reference to the description of the applied processing
      const std::vector<DataProcessing>& getDataProcessing() const;
      /// returns a mutable reference to the description of the applied processing
      std::vector<DataProcessing>& getDataProcessing();
      /// sets the description of the applied processing
      void setDataProcessing(const std::vector<DataProcessing>& data_processing);

    protected:
    	
    	SpectrumType type_;
    	String native_id_;
      String comment_;
      InstrumentSettings instrument_settings_;
      SourceFile source_file_;
      AcquisitionInfo acquisition_info_;
      std::vector<Precursor> precursors_;
      std::vector<Product> products_;
	    std::vector<PeptideIdentification> identification_;
		  std::vector<DataProcessing> data_processing_;
  };

	///Print the contents to a stream.
	OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const SpectrumSettings& spec);

} // namespace OpenMS

#endif // OPENMS_METADATA_SPECTRUMSETTINGS_H

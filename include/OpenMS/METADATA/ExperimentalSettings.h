// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_EXPERIMENTALSETTINGS_H
#define OPENMS_METADATA_EXPERIMENTALSETTINGS_H

#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/HPLC.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <vector>

namespace OpenMS 
{
	/**
		@brief Description of the experimental settings
		
		These settings are valid for the whole experiment. 
		See SpectrumSettings for settings which are spectrum specific.
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI ExperimentalSettings
  	: public MetaInfoInterface,
			public DocumentIdentifier
  {
    public:			

    	///Native ID type
    	enum NativeIDType
    	{
    		UNKNOWN,							///< Unknown native ID type
    		THERMO,								///< controller=xsd:nonNegativeInteger scan=xsd:positiveInteger
    		WATERS,								///< function=xsd:positiveInteger process=xsd:nonNegativeInteger scan=xsd:nonNegativeInteger
    		WIFF,									///< sample=xsd:nonNegativeInteger period=xsd:nonNegativeInteger cycle=xsd:nonNegativeInteger experiment=xsd:nonNegativeInteger
    		BRUKER_AGILENT,				///< scan=xsd:nonNegativeInteger
    		BRUKER_BAF,						///< scan=xsd:nonNegativeInteger
    		BRUKER_FID,						///< file=xsd:IDREF
    		MULTIPLE_PEAK_LISTS,	///< index=xsd:nonNegativeInteger @n Used for conversion of peak list files with multiple spectra, i.e. MGF, PKL, merged DTA files. Index is the spectrum number in the file, starting from 0.
    		SINGLE_PEAK_LIST,			///< file=xsd:IDREF @n The nativeID must be the same as the source file ID. Used for conversion of peak list files with one spectrum per file, typically folder of PKL or DTAs, each sourceFileRef is different.
    		SCAN_NUMBER,					///< scan=xsd:nonNegativeInteger @n Used for conversion from mzXML, or DTA folder where native scan numbers can be derived.
    		SPECTRUM_IDENTIFIER,	///< spectrum=xsd:nonNegativeInteger @n Used for conversion from mzData. The spectrum id attribute is referenced.
    		SIZE_OF_NATIVEIDTYPE
    	};
			/// Names of native ID types
			static const std::string NamesOfNativeIDType[SIZE_OF_NATIVEIDTYPE];

			///Constructor
      ExperimentalSettings();
      ///Copy constructor
      ExperimentalSettings(const ExperimentalSettings& source);
      ///Destructor
      ~ExperimentalSettings();
      
      ///Assignment operator
      ExperimentalSettings& operator= (const ExperimentalSettings& source);

      /// Equality operator
      bool operator== (const ExperimentalSettings& rhs) const;      
      /// Equality operator
      bool operator!= (const ExperimentalSettings& rhs) const;
			
			/// Returns the native ID type of the spectra
  		NativeIDType getNativeIDType() const;
			/// Sets the native ID type of the spectra
  		void setNativeIDType(NativeIDType type);

	    /// returns a const reference to the sample description
      const Sample& getSample() const;
      /// returns a mutable reference to the sample description
      Sample& getSample();
      /// sets the sample description
      void setSample(const Sample& sample);
			
			/// returns a const reference to the source date file
      const std::vector<SourceFile>& getSourceFiles() const;
      /// returns a mutable reference to the source date file
      std::vector<SourceFile>& getSourceFiles();
      /// sets the source date file
      void setSourceFiles(const std::vector<SourceFile>& source_files);
			
			/// returns a const reference to the list of contact persons
      const std::vector<ContactPerson>& getContacts() const;
      	/// returns a mutable reference to the list of contact persons
      std::vector<ContactPerson>& getContacts();
      	/// sets the list of contact persons
      void setContacts(const std::vector<ContactPerson>& contacts);
			
			/// returns a const reference to the MS instrument description
      const Instrument& getInstrument() const;
      /// returns a mutable reference to the MS instrument description
      Instrument& getInstrument();
      /// sets the MS instrument description
      void setInstrument(const Instrument& instrument);
			
			/// returns a const reference to the description of the applied processing 
      const std::vector<DataProcessing>& getDataProcessing() const;
      /// returns a mutable reference to the description of the applied processing 
      std::vector<DataProcessing>& getDataProcessing();
      /// sets the description of the applied processing 
      void setDataProcessing(const std::vector<DataProcessing>& processing_method);

			/// returns a const reference to the description of the HPLC run
      const HPLC& getHPLC() const;
      /// returns a mutable reference to the description of the HPLC run
      HPLC& getHPLC();
      /// sets the description of the HPLC run
      void setHPLC(const HPLC& hplc);
      
     	/// returns the date the experiment was performed
    	const DateTime& getDateTime() const;
    	/// sets the date the experiment was performed
      void setDateTime(const DateTime& date);   

			/// returns the free-text comment
      const String& getComment() const;
      /// sets the free-text comment
      void setComment(const String& comment);

		 	/// returns a const reference to the protein ProteinIdentification vector
		 	const std::vector<ProteinIdentification>& getProteinIdentifications() const;		 		    	
		 	/// returns a mutable reference to the protein ProteinIdentification vector
		  std::vector<ProteinIdentification>& getProteinIdentifications();		  
		  /// sets the protein ProteinIdentification vector
		  void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications);

    protected:
    	NativeIDType native_id_type_;
			Sample sample_;
			std::vector<SourceFile> source_files_;
			std::vector<ContactPerson> contacts_;
			Instrument instrument_;
		  std::vector<DataProcessing> data_processing_;
		  HPLC hplc_;
		  DateTime datetime_;
			String comment_;
			std::vector<ProteinIdentification> protein_identifications_;
  };

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const ExperimentalSettings& exp);

} // namespace OpenMS

#endif // OPENMS_METADATA_EXPERIMENTALSETTINGS_H

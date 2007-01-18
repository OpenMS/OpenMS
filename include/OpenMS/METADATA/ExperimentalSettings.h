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

#ifndef OPENMS_METADATA_EXPERIMENTALSETTINGS_H
#define OPENMS_METADATA_EXPERIMENTALSETTINGS_H

#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/METADATA/ProcessingMethod.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/HPLC.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <vector>

namespace OpenMS 
{
	/**
		@brief Description of the experimental settings
		
		These settings are valid for the whole experiment. 
		See SpectrumSettings for settings which are spectrum specific.
		
		@todo add several instances of software to capture processing in TOPP (Marc)
		
		@ingroup Metadata
	*/
  class ExperimentalSettings
  	: public MetaInfoInterface
  {
    public:
    	///Type of the experiment
    	enum ExperimentType 
    	{
    		UNKNOWN, ///< Unknown experiment type
    		MS, ///< MS experiment
    		MS_MS, ///< Tandem MS experiment
    		HPLC_MS, ///< HPLC-MS experiment
    		HPLC_MS_MS, ///< HPLC-MS experiment with tandem MS information
    		SIZE_OF_EXPERIMENTTYPE ///< Number of experiment types (can be used to iterate over the names)
    	};
    	///Names of experiment types                                 
    	static const std::string NamesOfExperimentType[SIZE_OF_EXPERIMENTTYPE];
			
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

	    /// returns a const reference to the sample description
      const Sample& getSample() const;
      /// returns a mutable reference to the sample description
      Sample& getSample();
      /// sets the sample description
      void setSample(const Sample& sample);
			
			/// returns a const reference to the source date file
      const SourceFile& getSourceFile() const;
      /// returns a mutable reference to the source date file
      SourceFile& getSourceFile();
      /// sets the source date file
      void setSourceFile(const SourceFile& source_file);
			
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
			
			/// returns a const reference to the software used for processing
      const Software& getSoftware() const;
      /// returns a mutable reference to the software used for processing
      Software& getSoftware();
      /// sets the software used for processing
      void setSoftware(const Software& software);
			
			/// returns a const reference to the description of the applied processing 
      const ProcessingMethod& getProcessingMethod() const;
      /// returns a mutable reference to the description of the applied processing 
      ProcessingMethod& getProcessingMethod();
      /// sets the description of the applied processing 
      void setProcessingMethod(const ProcessingMethod& processing_method);

			/// returns a const reference to the description of the HPLC run
      const HPLC& getHPLC() const;
      /// returns a mutable reference to the description of the HPLC run
      HPLC& getHPLC();
      /// sets the description of the HPLC run
      void setHPLC(const HPLC& hplc);

      /// returns the experiment type
    	ExperimentType getType() const;
    	/// sets the experiment type
      void setType(ExperimentType type);

     	/// returns the date the experiment was performed
    	const Date& getDate() const;
    	/// sets the date the experiment was performed
      void setDate(const Date& date);   

			/// returns the free-text comment
      const String& getComment() const;
      /// sets the free-text comment
      void setComment(const String& comment);

		 	/// returns a const reference to the protein identification vector
		 	const std::vector<ProteinIdentification>& getProteinIdentifications() const;		 		    	
		 	/// returns a mutable reference to the protein identification vector
		  std::vector<ProteinIdentification>& getProteinIdentifications();		  
		  /// sets the protein identification vector
		  void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications);
		  /// adds an identification to the identification vector
		  void addProteinIdentification(ProteinIdentification& protein_identification);

    protected:
			Sample sample_;
			SourceFile source_file_;
			std::vector<ContactPerson> contacts_;
			Instrument instrument_;
			Software software_;
		  ProcessingMethod processing_method_;
		  HPLC hplc_;
		  ExperimentType type_;
		  Date date_;
			String comment_;
			std::vector<ProteinIdentification> protein_identifications_;
  };

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const ExperimentalSettings& exp);

} // namespace OpenMS

#endif // OPENMS_METADATA_EXPERIMENTALSETTINGS_H

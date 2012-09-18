// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_CHROMATOGRAMSETTINGS_H
#define OPENMS_METADATA_CHROMATOGRAMSETTINGS_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Product.h>
#include <OpenMS/METADATA/DataProcessing.h>

#include <map>
#include <vector>

namespace OpenMS 
{
	/**
		@brief Representation of chromatogram settings, e.g. SRM/MRM chromatograms
		
		It contains the metadata about chromatogram specific instrument settings,
		acquisition settings, description of the meta values used in the peaks and precursor info.
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI ChromatogramSettings
  	: public MetaInfoInterface
  {
  	
    public:
   
			enum ChromatogramType
			{
				MASS_CHROMATOGRAM = 0,
				TOTAL_ION_CURRENT_CHROMATOGRAM,
				SELECTED_ION_CURRENT_CHROMATOGRAM,
				BASEPEAK_CHROMATOGRAM,
				SELECTED_ION_MONITORING_CHROMATOGRAM,
				SELECTED_REACTION_MONITORING_CHROMATOGRAM,
				ELECTROMAGNETIC_RADIATION_CHROMATOGRAM,
				ABSORPTION_CHROMATOGRAM,
				EMISSION_CHROMATOGRAM
			};

    	/// Constructor
      ChromatogramSettings();
      /// Copy constructor
      ChromatogramSettings(const ChromatogramSettings& source);
      /// Destructor
      virtual ~ChromatogramSettings();
      
      // Assignment operator
      ChromatogramSettings& operator= (const ChromatogramSettings& source);

      /// Equality operator
      bool operator== (const ChromatogramSettings& rhs) const;
      /// Equality operator
      bool operator!= (const ChromatogramSettings& rhs) const;


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
      const Precursor& getPrecursor() const;
      /// returns a mutable reference to the precursors
      Precursor& getPrecursor();
      /// sets the precursors
      void setPrecursor(const Precursor& precursor);

			/// returns a const reference to the products
      const Product& getProduct() const;
      /// returns a mutable reference to the products
      Product& getProduct();
      /// sets the products
      void setProduct(const Product& product);
			
			/// returns a const reference to the description of the applied processing
      const std::vector<DataProcessing>& getDataProcessing() const;
      /// returns a mutable reference to the description of the applied processing
      std::vector<DataProcessing>& getDataProcessing();
      /// sets the description of the applied processing
      void setDataProcessing(const std::vector<DataProcessing>& data_processing);

			/// returns the chromatogram type, e.g. a SRM chromatogram
			ChromatogramType getChromatogramType() const;

			/// sets the chromatogram type
			void setChromatogramType(ChromatogramType type);


    protected:
    	
    	String native_id_;
      String comment_;
      InstrumentSettings instrument_settings_;
      SourceFile source_file_;
      AcquisitionInfo acquisition_info_;
      Precursor precursor_;
      Product product_;
		  std::vector<DataProcessing> data_processing_;
			ChromatogramType type_;
  };

	///Print the contents to a stream.
	OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const ChromatogramSettings& spec);

} // namespace OpenMS

#endif // OPENMS_METADATA_CHROMATOGRAMSETTINGS_H

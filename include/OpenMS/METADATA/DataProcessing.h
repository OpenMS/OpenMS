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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_DATAPROCESSING_H
#define OPENMS_METADATA_DATAPROCESSING_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <set>

namespace OpenMS 
{
	/**
		@brief Descripton of the applied preprocessing steps
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI DataProcessing
  	: public MetaInfoInterface
  {
  	
    public:
    
    	//The different processing types
    	enum ProcessingAction
    	{
    		DATA_PROCESSING,          ///< General data processing (if no other term applies)
    		CHARGE_DECONVOLUTION,		  ///< Charge deconvolution
    		DEISOTOPING, 						  ///< Deisotoping
    		SMOOTHING, 							  ///< Smoothing of the signal to reduce noise
    		CHARGE_CALCULATION,       ///< Determination of the peak charge
    		PRECURSOR_RECALCULATION,	///< Recalculation of precursor m/z
    		BASELINE_REDUCTION, 		  ///< Baseline reduction
    		PEAK_PICKING, 					  ///< Peak picking (conversion from raw to peak data)
    		ALIGNMENT, 							  ///< Retention time alignment of different maps
    		CALIBRATION, 							///< Calibration of m/z positions
    		NORMALIZATION, 						///< Normalization of intensity values
    		FILTERING, 							  ///< Data filtering or extraction
    		QUANTITATION, 						///< Quantitation
    		FEATURE_GROUPING, 				///< %Feature grouping
    		IDENTIFICATION_MAPPING,		///< %Identification mapping
    		FORMAT_CONVERSION,        ///< General file format conversion (if no other term applies)
    		CONVERSION_MZDATA,			  ///< Convertion to mzData format
    		CONVERSION_MZML,				  ///< Conversion to mzML format
    		CONVERSION_MZXML,				  ///< Conversion to mzXML format
    		CONVERSION_DTA,           ///< Conversion to DTA format
    		SIZE_OF_PROCESSINGACTION
    	};
    	/// Names of inlet types
			static const std::string NamesOfProcessingAction[SIZE_OF_PROCESSINGACTION];
			
      /// Constructor
      DataProcessing();
      /// Copy construcor
      DataProcessing(const DataProcessing& source);
      /// Destructor
      ~DataProcessing();
      
      /// Assignement operator
      DataProcessing& operator= (const DataProcessing& source);

      /// Equality operator
      bool operator== (const DataProcessing& rhs) const;
      /// Equality operator
      bool operator!= (const DataProcessing& rhs) const;

			/// returns a const reference to the software used for processing
      const Software& getSoftware() const;
      /// returns a mutable reference to the software used for processing
      Software& getSoftware();
      /// sets the software used for processing
      void setSoftware(const Software& software);
      
    	/// returns a const reference to the applied processing actions
      const std::set<ProcessingAction>& getProcessingActions() const;
      /// returns a mutable reference to the description of the applied processing 
      std::set<ProcessingAction>& getProcessingActions();
      /// sets the description of the applied processing 
      void setProcessingActions(const std::set<ProcessingAction>& actions);

			/// returns the time of completition of the processing
	    const DateTime& getCompletionTime() const;
      /// sets the time of completition taking a DateTime object
      void setCompletionTime(const DateTime& completion_time);

    protected:

    	Software software_;
      std::set<ProcessingAction> processing_actions_;
      DateTime completion_time_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_DATAPROCESSING_H

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
// $Maintainer: $
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

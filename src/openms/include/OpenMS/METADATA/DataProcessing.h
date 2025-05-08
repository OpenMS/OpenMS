// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <set>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  /**
      @brief Description of the applied preprocessing steps

      @ingroup Metadata
  */
  class OPENMS_DLLAPI DataProcessing :
    public MetaInfoInterface
  {

public:

    //The different processing types
    enum ProcessingAction
    {
      DATA_PROCESSING,                ///< General data processing (if no other term applies)
      CHARGE_DECONVOLUTION,               ///< Charge deconvolution
      DEISOTOPING,                                ///< Deisotoping
      SMOOTHING,                                  ///< Smoothing of the signal to reduce noise
      CHARGE_CALCULATION,             ///< Determination of the peak charge
      PRECURSOR_RECALCULATION,          ///< Recalculation of precursor m/z
      BASELINE_REDUCTION,                 ///< Baseline reduction
      PEAK_PICKING,                           ///< Peak picking (conversion from raw to peak data)
      ALIGNMENT,                                  ///< Retention time alignment of different maps
      CALIBRATION,                                  ///< Calibration of m/z positions
      NORMALIZATION,                            ///< Normalization of intensity values
      FILTERING,                                  ///< Data filtering or extraction
      QUANTITATION,                             ///< Quantitation
      FEATURE_GROUPING,                     ///< %Feature grouping
      IDENTIFICATION_MAPPING,           ///< %Identification mapping
      FORMAT_CONVERSION,              ///< General file format conversion (if no other term applies)
      CONVERSION_MZDATA,                  ///< Conversion to mzData format
      CONVERSION_MZML,                        ///< Conversion to mzML format
      CONVERSION_MZXML,                       ///< Conversion to mzXML format
      CONVERSION_DTA,                 ///< Conversion to DTA format
      IDENTIFICATION,                 ///< Identification
      ION_MOBILITY_BINNING,           ///< Ion mobility binning (merging of spectra with similar IM values)
      SIZE_OF_PROCESSINGACTION
    };
    /// Names of inlet types
    static const std::string NamesOfProcessingAction[SIZE_OF_PROCESSINGACTION];

    /// Constructor
    DataProcessing() = default;
    /// Copy constructor
    DataProcessing(const DataProcessing&) = default;

    // note: we implement the move constructor ourselves due to a bug in MSVS
    // 2015/2017 which cannot produce a default move constructor for classes
    // that contain STL containers (other than vector).

    /// Move constructor
    DataProcessing(DataProcessing&&) noexcept;
    /// Destructor
    ~DataProcessing();

    /// Assignment operator
    DataProcessing& operator=(const DataProcessing&) = default;
    /// Move assignment operator
    DataProcessing& operator=(DataProcessing&&)& = default;

    /// Equality operator
    bool operator==(const DataProcessing& rhs) const;
    /// Equality operator
    bool operator!=(const DataProcessing& rhs) const;

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

    /// returns the time of completion of the processing
    const DateTime& getCompletionTime() const;
    /// sets the time of completion taking a DateTime object
    void setCompletionTime(const DateTime& completion_time);

protected:

    Software software_;
    std::set<ProcessingAction> processing_actions_;
    DateTime completion_time_;
  };

  typedef boost::shared_ptr<DataProcessing> DataProcessingPtr;
  typedef boost::shared_ptr<const DataProcessing> ConstDataProcessingPtr;

} // namespace OpenMS

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

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
  class OPENMS_DLLAPI ChromatogramSettings :
    public MetaInfoInterface
  {

public:

    /// List of chromatogram names, e.g., SELECTED_REACTION_MONITORING_CHROMATOGRAM.
    /// Actual names can be accessed using the ChromatogramNames[] array
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
      EMISSION_CHROMATOGRAM,
      SIZE_OF_CHROMATOGRAM_TYPE // last entry!
    };

    /// Names of chromatogram types corresponding to enum ChromatogramType
    static const char * const ChromatogramNames[SIZE_OF_CHROMATOGRAM_TYPE+1]; // avoid string[], since it gets copied onto heap on initialization

    /// Constructor
    ChromatogramSettings();
    /// Copy constructor
    ChromatogramSettings(const ChromatogramSettings &) = default;
    /// Move constructor
    ChromatogramSettings(ChromatogramSettings&&) = default;
    /// Destructor
    virtual ~ChromatogramSettings();

    // Assignment operator
    ChromatogramSettings & operator=(const ChromatogramSettings &) = default;
    /// Move assignment operator
    ChromatogramSettings& operator=(ChromatogramSettings&&) & = default;

    /// Equality operator
    bool operator==(const ChromatogramSettings & rhs) const;
    /// Equality operator
    bool operator!=(const ChromatogramSettings & rhs) const;

    /// returns the native identifier for the spectrum, used by the acquisition software.
    const String & getNativeID() const;
    /// sets the native identifier for the spectrum, used by the acquisition software.
    void setNativeID(const String & native_id);

    /// returns the free-text comment
    const String & getComment() const;
    /// sets the free-text comment
    void setComment(const String & comment);

    /// returns a const reference to the instrument settings of the current spectrum
    const InstrumentSettings & getInstrumentSettings() const;
    /// returns a mutable reference to the instrument settings of the current spectrum
    InstrumentSettings & getInstrumentSettings();
    /// sets the instrument settings of the current spectrum
    void setInstrumentSettings(const InstrumentSettings & instrument_settings);

    /// returns a const reference to the acquisition info
    const AcquisitionInfo & getAcquisitionInfo() const;
    /// returns a mutable reference to the acquisition info
    AcquisitionInfo & getAcquisitionInfo();
    /// sets the acquisition info
    void setAcquisitionInfo(const AcquisitionInfo & acquisition_info);

    /// returns a const reference to the source file
    const SourceFile & getSourceFile() const;
    /// returns a mutable reference to the source file
    SourceFile & getSourceFile();
    /// sets the source file
    void setSourceFile(const SourceFile & source_file);

    /// returns a const reference to the precursors
    const Precursor & getPrecursor() const;
    /// returns a mutable reference to the precursors
    Precursor & getPrecursor();
    /// sets the precursors
    void setPrecursor(const Precursor & precursor);

    /// returns a const reference to the products
    const Product & getProduct() const;
    /// returns a mutable reference to the products
    Product & getProduct();
    /// sets the products
    void setProduct(const Product & product);

    /// returns the chromatogram type, e.g. a SRM chromatogram
    ChromatogramType getChromatogramType() const;

    /// sets the chromatogram type
    void setChromatogramType(ChromatogramType type);

    /// sets the description of the applied processing
    void setDataProcessing(const std::vector< DataProcessingPtr > & data_processing);

    /// returns a mutable reference to the description of the applied processing
    std::vector< DataProcessingPtr > & getDataProcessing();

    /// returns a const reference to the description of the applied processing
    const std::vector< boost::shared_ptr<const DataProcessing > > getDataProcessing() const;

protected:

    String native_id_;
    String comment_;
    InstrumentSettings instrument_settings_;
    SourceFile source_file_;
    AcquisitionInfo acquisition_info_;
    Precursor precursor_;
    Product product_;
    std::vector< DataProcessingPtr > data_processing_;
    ChromatogramType type_;
  };

  ///Print the contents to a stream.
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const ChromatogramSettings & spec);

} // namespace OpenMS


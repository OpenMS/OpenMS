// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

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
  class OPENMS_DLLAPI SpectrumSettings :
    public MetaInfoInterface
  {

public:

    ///Spectrum peak type
    enum SpectrumType
    {
      UNKNOWN,          ///< Unknown spectrum type
      CENTROID,         ///< centroid data or stick data
      PROFILE,          ///< profile data
      SIZE_OF_SPECTRUMTYPE
    };
    /// Names of spectrum types
    static const std::string NamesOfSpectrumType[SIZE_OF_SPECTRUMTYPE];

    /// Constructor
    SpectrumSettings();
    /// Copy constructor
    SpectrumSettings(const SpectrumSettings &) = default;
    /// Move constructor
    SpectrumSettings(SpectrumSettings&&) = default;
    /// Destructor
    ~SpectrumSettings();

    // Assignment operator
    SpectrumSettings & operator=(const SpectrumSettings &) = default;
    /// Move assignment operator
    SpectrumSettings& operator=(SpectrumSettings&&) & = default;

    /// Equality operator
    bool operator==(const SpectrumSettings & rhs) const;
    /// Equality operator
    bool operator!=(const SpectrumSettings & rhs) const;

    /// merge another spectrum setting into this one (data is usually appended, except for spectrum type which needs to be unambiguous to be kept)
    void unify(const SpectrumSettings & rhs);

    ///returns the spectrum type (centroided (PEAKS) or profile data (RAW))
    SpectrumType getType() const;
    ///sets the spectrum type
    void setType(SpectrumType type);

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
    const std::vector<Precursor> & getPrecursors() const;
    /// returns a mutable reference to the precursors
    std::vector<Precursor> & getPrecursors();
    /// sets the precursors
    void setPrecursors(const std::vector<Precursor> & precursors);

    /// returns a const reference to the products
    const std::vector<Product> & getProducts() const;
    /// returns a mutable reference to the products
    std::vector<Product> & getProducts();
    /// sets the products
    void setProducts(const std::vector<Product> & products);

    /// returns a const reference to the PeptideIdentification vector
    const std::vector<PeptideIdentification> & getPeptideIdentifications() const;
    /// returns a mutable reference to the PeptideIdentification vector
    std::vector<PeptideIdentification> & getPeptideIdentifications();
    /// sets the PeptideIdentification vector
    void setPeptideIdentifications(const std::vector<PeptideIdentification> & identifications);

    /// sets the description of the applied processing
    void setDataProcessing(const std::vector< DataProcessingPtr > & data_processing);

    /// returns a mutable reference to the description of the applied processing
    std::vector< DataProcessingPtr > & getDataProcessing();

    /// returns a const reference to the description of the applied processing
    const std::vector< boost::shared_ptr<const DataProcessing > > getDataProcessing() const;

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
    std::vector< DataProcessingPtr > data_processing_;
  };

  ///Print the contents to a stream.
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const SpectrumSettings & spec);

} // namespace OpenMS


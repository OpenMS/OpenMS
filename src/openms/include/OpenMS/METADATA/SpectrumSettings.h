// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Product.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/HashUtils.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>

#include <functional>
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
    enum class SpectrumType
    {
      UNKNOWN,          ///< Unknown spectrum type
      CENTROID,         ///< centroid data or stick data
      PROFILE,          ///< profile data
      SIZE_OF_SPECTRUMTYPE
    };
    /// Names of spectrum types
    static const std::string NamesOfSpectrumType[static_cast<size_t>(SpectrumType::SIZE_OF_SPECTRUMTYPE)];

    /// returns all spectrum type names known to OpenMS
    static StringList getAllNamesOfSpectrumType();

    /// Convert a SpectrumType enum to String
    /// @throws Exception::InvalidValue if @p type is SIZE_OF_SPECTRUMTYPE
    static const std::string& spectrumTypeToString(SpectrumType type);

    /// Convert an entry in NamesOfSpectrumType[] to SpectrumType enum
    /// @throws Exception::InvalidValue if @p name is not contained in NamesOfSpectrumType[]
    static SpectrumType toSpectrumType(const std::string& name);

    /// Constructor
    SpectrumSettings() = default;
    /// Copy constructor
    SpectrumSettings(const SpectrumSettings &) = default;
    /// Move constructor
    SpectrumSettings(SpectrumSettings&&) noexcept = default;
    /// Destructor
    ~SpectrumSettings() noexcept = default;

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

    /// @brief sets the IMFormat of the spectrum
    /// @param[in] im_type
    void setIMFormat(const IMFormat& im_type);

    /// @brief returns the IMFormat of the spectrum if set. Otherwise UNKNOWN (default).
    ///
    /// Note: If UNKNOWN, use IMFormat::determineIMFormat to determine the IMFormat based on the data.
    /// @return IMFormat of the spectrum
    IMFormat getIMFormat() const;

    /// @brief sets the IM peak type (profile vs centroided in the IM dimension)
    /// @param[in] im_peak_type the IM processing state to set
    void setIMPeakType(IMPeakType im_peak_type);

    /// @brief returns the IM peak type of the spectrum if set. Otherwise UNKNOWN (default).
    /// @return IMPeakType of the spectrum
    IMPeakType getIMPeakType() const;

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

    /// sets the description of the applied processing
    void setDataProcessing(const std::vector< DataProcessingPtr > & data_processing);

    /// returns a mutable reference to the description of the applied processing
    std::vector< DataProcessingPtr > & getDataProcessing();

    /// returns a const reference to the description of the applied processing
    const std::vector< std::shared_ptr<const DataProcessing > > getDataProcessing() const;

protected:

    SpectrumType type_ = SpectrumType::UNKNOWN;
    IMFormat im_type_ = IMFormat::UNKNOWN;
    IMPeakType im_peak_type_ = IMPeakType::UNKNOWN;
    String native_id_;
    String comment_;
    InstrumentSettings instrument_settings_;
    SourceFile source_file_;
    AcquisitionInfo acquisition_info_;
    std::vector<Precursor> precursors_;
    std::vector<Product> products_;
    std::vector< DataProcessingPtr > data_processing_;
  };

  ///Print the contents to a stream.
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const SpectrumSettings & spec);

} // namespace OpenMS

namespace std
{
  /**
   * @brief Hash function for OpenMS::SpectrumSettings.
   *
   * Hashes all fields used in operator==:
   * - MetaInfoInterface (via getKeys/getMetaValue)
   * - type_ (SpectrumType enum)
   * - native_id_ (String)
   * - comment_ (String)
   * - instrument_settings_ (InstrumentSettings - hashes scan_mode, zoom_scan, polarity, scan_windows)
   * - acquisition_info_ (AcquisitionInfo - hashes method_of_combination and acquisitions)
   * - source_file_ (SourceFile - hashes name, path, file_size, file_type, checksum, native_id_type)
   * - precursors_ (vector<Precursor> - hashes each precursor's key fields)
   * - products_ (vector<Product> - uses std::hash<Product>)
   * - data_processing_ (vector<DataProcessingPtr> - hashes dereferenced contents)
   *
   * @note This hash is consistent with operator==. Since operator== compares
   * all fields including MetaInfoInterface and complex nested types, this hash
   * includes representative fields from each. For nested types without their own
   * std::hash, we hash their primary identifying fields.
   */
  template<>
  struct hash<OpenMS::SpectrumSettings>
  {
    std::size_t operator()(const OpenMS::SpectrumSettings& s) const noexcept
    {
      std::size_t seed = 0;

      // Hash type_ (enum)
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(s.getType())));

      // Hash native_id_ (String)
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(s.getNativeID()));

      // Hash comment_ (String)
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(s.getComment()));

      // Hash instrument_settings_ (all fields from operator==)
      const auto& is = s.getInstrumentSettings();
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(is.getScanMode())));
      OpenMS::hash_combine(seed, OpenMS::hash_int(is.getZoomScan() ? 1 : 0));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(is.getPolarity())));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(is.getScanWindows().size())));
      for (const auto& sw : is.getScanWindows())
      {
        OpenMS::hash_combine(seed, OpenMS::hash_float(sw.begin));
        OpenMS::hash_combine(seed, OpenMS::hash_float(sw.end));
      }
      OpenMS::hash_combine(seed, std::hash<OpenMS::MetaInfoInterface>{}(is));

      // Hash acquisition_info_ (all fields from operator==)
      const auto& ai = s.getAcquisitionInfo();
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(ai.getMethodOfCombination()));
      OpenMS::hash_combine(seed, std::hash<OpenMS::MetaInfoInterface>{}(ai));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(ai.size())));
      for (const auto& acq : ai)
      {
        OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(acq.getIdentifier()));
        OpenMS::hash_combine(seed, std::hash<OpenMS::MetaInfoInterface>{}(acq));
      }

      // Hash source_file_ (all fields from operator==, including CVTermList base)
      const auto& sf = s.getSourceFile();
      OpenMS::hash_combine(seed, OpenMS::hashCVTermList(sf));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getNameOfFile()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getPathToFile()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(sf.getFileSize()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getFileType()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getChecksum()));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(sf.getChecksumType())));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getNativeIDType()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getNativeIDTypeAccession()));

      // Hash precursors_ (using std::hash<Precursor>)
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(s.getPrecursors().size())));
      for (const auto& p : s.getPrecursors())
      {
        OpenMS::hash_combine(seed, std::hash<OpenMS::Precursor>{}(p));
      }

      // Hash products_ (using std::hash<Product>)
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(s.getProducts().size())));
      for (const auto& p : s.getProducts())
      {
        OpenMS::hash_combine(seed, std::hash<OpenMS::Product>{}(p));
      }

      // Hash data_processing_ (dereferenced contents)
      const auto dp = s.getDataProcessing();
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(dp.size())));
      for (const auto& dp_ptr : dp)
      {
        if (dp_ptr)
        {
          // Hash Software name and version
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(dp_ptr->getSoftware().getName()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(dp_ptr->getSoftware().getVersion()));
          // Hash processing actions
          for (const auto& action : dp_ptr->getProcessingActions())
          {
            OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(action)));
          }
          // Hash completion time
          OpenMS::hash_combine(seed, std::hash<OpenMS::DateTime>{}(dp_ptr->getCompletionTime()));
        }
      }

      // Hash MetaInfoInterface base class
      OpenMS::hash_combine(seed, std::hash<OpenMS::MetaInfoInterface>{}(s));

      return seed;
    }
  };
} // namespace std


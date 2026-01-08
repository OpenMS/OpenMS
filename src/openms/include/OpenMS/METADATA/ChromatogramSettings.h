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
#include <OpenMS/CONCEPT/HashUtils.h>

#include <functional>
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
    const std::vector< std::shared_ptr<const DataProcessing > > getDataProcessing() const;

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

// Hash function specialization for ChromatogramSettings
namespace std
{
  /**
   * @brief Hash function for OpenMS::ChromatogramSettings.
   *
   * Hashes all fields used in operator==:
   * - MetaInfoInterface (base class - meta values)
   * - native_id_ (string)
   * - comment_ (string)
   * - instrument_settings_ (scan mode, zoom scan, polarity, scan windows, and MetaInfo)
   * - acquisition_info_ (method of combination, acquisitions, and MetaInfo)
   * - source_file_ (file properties and CVTermList)
   * - precursor_ (activation methods, windows, drift time, charge, and Peak1D/CVTermList)
   * - product_ (uses std::hash<Product>)
   * - data_processing_ (software, actions, completion time for each element)
   * - type_ (chromatogram type enum)
   *
   * @note This hash implementation hashes the key identifying fields of nested types
   *       and includes MetaInfoInterface. Consistent with operator==.
   */
  template<>
  struct hash<OpenMS::ChromatogramSettings>
  {
    std::size_t operator()(const OpenMS::ChromatogramSettings& cs) const noexcept
    {
      std::size_t seed = 0;

      // Hash native_id_ and comment_
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cs.getNativeID()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cs.getComment()));

      // Hash instrument_settings_
      const auto& is = cs.getInstrumentSettings();
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(is.getScanMode())));
      OpenMS::hash_combine(seed, OpenMS::hash_int(is.getZoomScan() ? 1 : 0));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(is.getPolarity())));
      // Hash scan windows
      const auto& scan_windows = is.getScanWindows();
      OpenMS::hash_combine(seed, OpenMS::hash_int(scan_windows.size()));
      for (const auto& sw : scan_windows)
      {
        OpenMS::hash_combine(seed, OpenMS::hash_float(sw.begin));
        OpenMS::hash_combine(seed, OpenMS::hash_float(sw.end));
      }

      // Hash acquisition_info_
      const auto& ai = cs.getAcquisitionInfo();
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(ai.getMethodOfCombination()));
      OpenMS::hash_combine(seed, OpenMS::hash_int(ai.size()));
      for (const auto& acq : ai)
      {
        OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(acq.getIdentifier()));
      }

      // Hash source_file_
      const auto& sf = cs.getSourceFile();
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getNameOfFile()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getPathToFile()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(static_cast<double>(sf.getFileSize())));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getFileType()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getChecksum()));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(sf.getChecksumType())));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getNativeIDType()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(sf.getNativeIDTypeAccession()));

      // Hash precursor_
      const auto& prec = cs.getPrecursor();
      // Hash activation methods
      const auto& am = prec.getActivationMethods();
      OpenMS::hash_combine(seed, OpenMS::hash_int(am.size()));
      for (const auto& method : am)
      {
        OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(method)));
      }
      OpenMS::hash_combine(seed, OpenMS::hash_float(prec.getActivationEnergy()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(prec.getIsolationWindowLowerOffset()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(prec.getIsolationWindowUpperOffset()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(prec.getDriftTime()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(prec.getDriftTimeWindowLowerOffset()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(prec.getDriftTimeWindowUpperOffset()));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(prec.getDriftTimeUnit())));
      OpenMS::hash_combine(seed, OpenMS::hash_int(prec.getCharge()));
      // Hash possible charge states
      const auto& pcs = prec.getPossibleChargeStates();
      OpenMS::hash_combine(seed, OpenMS::hash_int(pcs.size()));
      for (const auto& charge : pcs)
      {
        OpenMS::hash_combine(seed, OpenMS::hash_int(charge));
      }
      // Hash Peak1D part (MZ and intensity)
      OpenMS::hash_combine(seed, OpenMS::hash_float(prec.getMZ()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(prec.getIntensity()));

      // Hash product_ using existing std::hash<Product>
      OpenMS::hash_combine(seed, std::hash<OpenMS::Product>{}(cs.getProduct()));

      // Hash data_processing_
      const auto dp = cs.getDataProcessing();
      OpenMS::hash_combine(seed, OpenMS::hash_int(dp.size()));
      for (const auto& proc_ptr : dp)
      {
        if (proc_ptr)
        {
          // Hash software using std::hash<Software>
          OpenMS::hash_combine(seed, std::hash<OpenMS::Software>{}(proc_ptr->getSoftware()));
          // Hash processing actions
          const auto& actions = proc_ptr->getProcessingActions();
          OpenMS::hash_combine(seed, OpenMS::hash_int(actions.size()));
          for (const auto& action : actions)
          {
            OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(action)));
          }
          // Hash completion time using std::hash<DateTime>
          OpenMS::hash_combine(seed, std::hash<OpenMS::DateTime>{}(proc_ptr->getCompletionTime()));
        }
      }

      // Hash type_
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(cs.getChromatogramType())));

      // Hash MetaInfoInterface base class
      OpenMS::hash_combine(seed, std::hash<OpenMS::MetaInfoInterface>{}(cs));

      return seed;
    }
  };
} // namespace std


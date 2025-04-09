// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config

#include <OpenMS/PROCESSING/MISC/DataFilters.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/VISUAL/INTERFACES/IPeptideIds.h>


namespace OpenMS
{

  /**
      @brief Base class to store either the currently visible or all data of a canvas
  */
  class OPENMS_GUI_DLLAPI LayerStoreData
  {
  public:
    LayerStoreData(FileTypeList supported_storage_formats) :
        storage_formats_(supported_storage_formats)
    {
    }

    /// virtual D'tor for proper cleanup of derived classes' members
    virtual ~LayerStoreData() = default;

    /// Which formats are supported when writing the file?
    FileTypeList getSupportedFileFormats() const
    {
      return storage_formats_;
    }

    /// Save the internal data to a file. The @p filename's suffix determines the file format. It must be one of getSupportedFileFormats() or UNKNOWN.
    /// If the @p filename's suffix is unknown, the first item from getSupportedFileFormats() determines the storage format.
    /// @param filename A relative or absolute path+filename. Its suffix determines the format.
    /// @param lt Show a progress bar in the GUI?
    /// @throw Exception::UnableToCreateFile if the extension of @p filename is neither in getSupportedFileFormats() nor UNKNOWN.
    virtual void saveToFile(const String& filename, const ProgressLogger::LogType lt) const = 0;

  protected:

    /// extracts the supported extension (converting UNKNOWN to first item in storage_formats_) or throws an Exception::UnableToCreateFile
    FileTypes::Type getSupportedExtension_(const String& filename) const;

    FileTypeList storage_formats_; ///< file formats which can hold the data from the layer; The first item should be the preferred/default format
  };

  /**
    @brief Visitor which can save a visible piece of data; subsequently the data can be stored to a file.
  */
  class OPENMS_GUI_DLLAPI LayerStoreDataPeakMapVisible
   : public LayerStoreData
  {
  public:
    LayerStoreDataPeakMapVisible() :
      LayerStoreData(FileTypeList({FileTypes::MZML, FileTypes::MZDATA, FileTypes::MZXML}))
    {}

    // docu in base class
    void saveToFile(const String& filename, const ProgressLogger::LogType lt) const override;

    /**
     * \brief Stores data from a 1D canvas and remembers the data internally
     * \param spec The spectrum to store
     * \param visible_range Restricts m/z (and intensity)
     * \param layer_filters Remove all peaks not passing this filter
     */
    void storeVisibleSpectrum(const MSSpectrum& spec, const RangeAllType& visible_range, const DataFilters& layer_filters);
    /**
     * \brief Stores data from a 1D canvas and remembers the data internally
     * \param chrom The chromatogram to store
     * \param visible_range Restricts RT (and intensity)
     * \param layer_filters Remove all peaks not passing this filter
     */
    void storeVisibleChromatogram(const MSChromatogram& chrom, const RangeAllType& visible_range, const DataFilters& layer_filters);

    /// analog to storeVisibleSpectrum()
    void storeVisibleExperiment(const PeakMap& exp, const RangeAllType& visible_range, const DataFilters& layer_filters);

  private:
    PeakMap pm_; ///< the filtered data; used when saveToFile() is called
  };

  /**
    @brief Visitor which can save a full experiment; subsequently the data can be stored to a file.

    Since only a pointer is stored internally, make sure the lifetime of the PeakMap exceeds this visitor.
  */
  class OPENMS_GUI_DLLAPI LayerStoreDataPeakMapAll : public LayerStoreData
  {
  public:
    LayerStoreDataPeakMapAll() : LayerStoreData(FileTypeList({FileTypes::MZML, FileTypes::MZDATA, FileTypes::MZXML}))
    {
    }

    // docu in base class
    void saveToFile(const String& filename, const ProgressLogger::LogType lt) const override;

    void storeFullExperiment(const PeakMap& exp);

  private:
    const PeakMap* full_exp_ = nullptr; ///< pointer to the full data, when storeFullExperiment() was called
  };

  /**
    @brief Visitor which can save a visible piece of data; subsequently the data can be stored to a file.
  */
  class OPENMS_GUI_DLLAPI LayerStoreDataFeatureMapVisible : public LayerStoreData
  {
  public:
    LayerStoreDataFeatureMapVisible() : LayerStoreData(FileTypeList({FileTypes::FEATUREXML}))
    {
    }

    // docu in base class
    void saveToFile(const String& filename, const ProgressLogger::LogType lt) const override;

    void storeVisibleFM(const FeatureMap& fm, const RangeAllType& visible_range, const DataFilters& layer_filters);

  private:
    FeatureMap fm_; ///< the filtered data; used when saveToFile() is called
  };

 /**
  @brief Visitor which can save a full FeatureMap; subsequently the data can be stored to a file.

  Since only a pointer is stored internally, make sure the lifetime of the FeatureMap exceeds this visitor.
*/
  class OPENMS_GUI_DLLAPI LayerStoreDataFeatureMapAll : public LayerStoreData
  {
  public:
    LayerStoreDataFeatureMapAll() : LayerStoreData(FileTypeList({FileTypes::FEATUREXML}))
    {
    }

    // docu in base class
    void saveToFile(const String& filename, const ProgressLogger::LogType lt) const override;

    void storeFullFM(const FeatureMap& fm);

  private:
    const FeatureMap* full_fm_ = nullptr; ///< pointer to the full data, when storeFullExperiment() was called
  };

  /**
    @brief Visitor which can save a visible piece of data; subsequently the data can be stored to a file.
  */
  class OPENMS_GUI_DLLAPI LayerStoreDataConsensusMapVisible : public LayerStoreData
  {
  public:
    LayerStoreDataConsensusMapVisible() : LayerStoreData(FileTypeList({FileTypes::CONSENSUSXML}))
    {
    }

    // docu in base class
    void saveToFile(const String& filename, const ProgressLogger::LogType lt) const override;

    void storeVisibleCM(const ConsensusMap& cm, const RangeAllType& visible_range, const DataFilters& layer_filters);

  private:
    ConsensusMap cm_; ///< the filtered data; used when saveToFile() is called
  };

  /**
   @brief Visitor which can save a full ConsensusMap; subsequently the data can be stored to a file.

   Since only a pointer is stored internally, make sure the lifetime of the ConsensusMap exceeds this visitor.
 */
  class OPENMS_GUI_DLLAPI LayerStoreDataConsensusMapAll : public LayerStoreData
  {
  public:
    LayerStoreDataConsensusMapAll() : LayerStoreData(FileTypeList({FileTypes::CONSENSUSXML}))
    {
    }

    // docu in base class
    void saveToFile(const String& filename, const ProgressLogger::LogType lt) const override;

    void storeFullCM(const ConsensusMap& cm);

  private:
    const ConsensusMap* full_cm_ = nullptr; ///< pointer to the full data, when storeFullExperiment() was called
  };



  /**
    @brief Visitor which can save a visible piece of data; subsequently the data can be stored to a file.
  */
  class OPENMS_GUI_DLLAPI LayerStoreDataIdentVisible : public LayerStoreData
  {
  public:
    LayerStoreDataIdentVisible() : LayerStoreData(FileTypeList({FileTypes::IDXML}))
    {
    }

    // docu in base class
    void saveToFile(const String& filename, const ProgressLogger::LogType lt) const override;

    void storeVisibleIdent(const IPeptideIds::PepIds& ids, const RangeAllType& visible_range, const DataFilters& layer_filters);

  private:
    IPeptideIds::PepIds ids_; ///< the filtered data; used when saveToFile() is called
  };

  /**
   @brief Visitor which can save a full set of Identifications; subsequently the data can be stored to a file.

   Since only a pointer is stored internally, make sure the lifetime of the Identifications exceeds this visitor.
 */
  class OPENMS_GUI_DLLAPI LayerStoreDataIdentAll : public LayerStoreData
  {
  public:
    LayerStoreDataIdentAll() : LayerStoreData(FileTypeList({FileTypes::IDXML}))
    {
    }

    // docu in base class
    void saveToFile(const String& filename, const ProgressLogger::LogType lt) const override;

    void storeFullIdent(const IPeptideIds::PepIds& ids);

  private:
    const IPeptideIds::PepIds* full_ids_ = nullptr; ///< pointer to the full data, when storeFullExperiment() was called
  };

} // namespace OpenMS
// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/VISUAL/LogWindow.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>
#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>

#include <boost/shared_ptr.hpp>

#include <vector>
#include <bitset>

class QWidget;

namespace OpenMS
{

  class OnDiscMSExperiment;
  class OSWData;

  /**
  @brief Class that stores the data for one layer

  The data for a layer can be peak data, feature data (feature, consensus),
  chromatogram or peptide identification data. 

  For 2D and 3D data, the data is generally accessible through getPeakData()
  while features are accessible through getFeatureMap() and getConsensusMap().
  For 1D data, the current spectrum must be accessed through
  getCurrentSpectrum().

  Peak data is stored using a shared pointer to an MSExperiment data structure
  as well as a shared pointer to a OnDiscMSExperiment data structure. Note that
  the actual data may not be in memory as this is not efficient for large files
  and therefore may have to be retrieved from disk on-demand. 

  @note The spectrum for 1D viewing retrieved through getCurrentSpectrum() is a
  copy of the actual raw data and *different* from the one retrieved through
  getPeakData()[index]. Any changes to applied to getCurrentSpectrum() are
  non-persistent and will be gone the next time the cache is updated.
  Persistent changes can be applied to getPeakDataMuteable() and will be
  available on the next cache update.

  @note Layer is mainly used as a member variable of PlotCanvas which holds
  a vector of LayerData objects.

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerData
  {
public:
    /** @name Type definitions */
    //@{
    /// Dataset types.
    /// Order in the enum determines the order in which layer types are drawn.
    enum DataType
    {
      DT_PEAK,            ///< Spectrum profile or centroided data
      DT_CHROMATOGRAM,    ///< Chromatogram data
      DT_FEATURE,         ///< Feature data
      DT_CONSENSUS,       ///< Consensus feature data
      DT_IDENT,           ///< Peptide identification data
      DT_UNKNOWN          ///< Undefined data type indicating an error
    };

    /// Flags that determine which information is shown.
    enum Flags
    {
      F_HULL,          ///< Features: Overall convex hull
      F_HULLS,         ///< Features: Convex hulls of single mass traces
      F_UNASSIGNED,    ///< Features: Unassigned peptide hits
      P_PRECURSORS,    ///< Peaks: Mark precursor peaks of MS/MS scans
      P_PROJECTIONS,   ///< Peaks: Show projections
      C_ELEMENTS,      ///< Consensus features: Show elements
      I_PEPTIDEMZ,     ///< Identifications: m/z source
      I_LABELS,        ///< Identifications: Show labels (not sequences)
      SIZE_OF_FLAGS
    };

    /// Actual state of each flag
    std::bitset<SIZE_OF_FLAGS> flags;

    /// Label used in visualization
    enum LabelType
    {
      L_NONE,          ///< No label is displayed
      L_INDEX,         ///< The element number is used
      L_META_LABEL,    ///< The 'label' meta information is used
      L_ID,            ///< The best peptide hit of the first identification run is used
      L_ID_ALL,        ///< All peptide hits of the first identification run are used
      SIZE_OF_LABEL_TYPE
    };

    /// Label names
    static const std::string NamesOfLabelType[SIZE_OF_LABEL_TYPE];

    /// Features
    typedef FeatureMap FeatureMapType;

    /// SharedPtr on feature map
    typedef boost::shared_ptr<FeatureMap > FeatureMapSharedPtrType;

    /// consensus features
    typedef ConsensusMap ConsensusMapType;

    /// SharedPtr on consensus features
    typedef boost::shared_ptr<ConsensusMap> ConsensusMapSharedPtrType;

    /// Main data type (experiment)
    typedef PeakMap ExperimentType;

    /// SharedPtr on MSExperiment
    typedef boost::shared_ptr<ExperimentType> ExperimentSharedPtrType;

    typedef boost::shared_ptr<const ExperimentType> ConstExperimentSharedPtrType;

    /// SharedPtr on On-Disc MSExperiment
    typedef boost::shared_ptr<OnDiscMSExperiment> ODExperimentSharedPtrType;

    /// SharedPtr on OSWData
    typedef boost::shared_ptr<OSWData> OSWDataSharedPtrType;

    //@}

    /// Default constructor
    LayerData();

    /// no Copy-ctor (should not be needed)
    LayerData(const LayerData& ld) = delete;
    /// no assignment operator (should not be needed)
    LayerData& operator=(const LayerData& ld) = delete;

    /// move Ctor
    LayerData(LayerData&& ld) = default;

    /// move assignment
    LayerData& operator=(LayerData&& ld) = default;

    /// Returns a const reference to the current feature data
    const FeatureMapSharedPtrType & getFeatureMap() const
    {
      return features_;
    }

    /// Returns a const reference to the current feature data
    FeatureMapSharedPtrType & getFeatureMap()
    {
      return features_;
    }

    /// Returns a const reference to the consensus feature data
    const ConsensusMapSharedPtrType & getConsensusMap() const
    {
      return consensus_map_;
    }

    /// Returns current consensus map (mutable)
    ConsensusMapSharedPtrType & getConsensusMap()
    {
      return consensus_map_;
    }

    /**
    @brief Returns a const reference to the current in-memory peak data

    @note Depending on the caching strategy (on-disk or in-memory), all or some
    spectra may have zero size and contain only meta data since peak data is
    cached on disk.

    @note Do *not* use this function to access the current spectrum for the 1D view, use getCurrentSpectrum() instead.
    */
    const ConstExperimentSharedPtrType getPeakData() const;

    /**
    @brief Returns a mutable reference to the current in-memory peak data

    @note Depending on the caching strategy (on-disk or in-memory), all or some
    spectra may have zero size and contain only meta data since peak data is
    cached on disk.

    @note Do *not* use this function to access the current spectrum for the 1D view, use getCurrentSpectrum() instead.
    */
    const ExperimentSharedPtrType & getPeakDataMuteable() {return peak_map_;}

    /**
    @brief Set the current in-memory peak data
    */
    void setPeakData(ExperimentSharedPtrType p)
    {
      peak_map_ = p;
      updateCache_();
    }

    /// Set the current on-disc data
    void setOnDiscPeakData(ODExperimentSharedPtrType p)
    {
      on_disc_peaks = p;
    }

    /// Returns a mutable reference to the on-disc data
    const ODExperimentSharedPtrType & getOnDiscPeakData() const
    {
      return on_disc_peaks;
    }

    /// Returns a mutable reference to the current chromatogram data
    const ExperimentSharedPtrType & getChromatogramData() const
    {
      return chromatogram_map_;
    }

    /// Returns a mutable reference to the current chromatogram data
    ExperimentSharedPtrType & getChromatogramData()
    {
      return chromatogram_map_;
    }

    /// get annotation (e.g. to build a hierachical ID View)
    /// Not const, because we might have incomplete data, which needs to be loaded from sql source
    OSWDataSharedPtrType& getChromatogramAnnotation();

    /// get annotation (e.g. to build a hierachical ID View)
    /// Not actually const (only the pointer, not the data), because we might have incomplete data, which needs to be loaded from sql source
    const OSWDataSharedPtrType& getChromatogramAnnotation() const;

    /// add annotation from an OSW sqlite file.
    void setChromatogramAnnotation(OSWData&& data);

    /// add peptide identifications to the layer
    /// Only supported for DT_PEAK, DT_FEATURE and DT_CONSENSUS.
    /// Will return false otherwise.
    bool annotate(const std::vector<PeptideIdentification>& identifications,
                  const std::vector<ProteinIdentification>& protein_identifications);


    /// Returns a const reference to the annotations of the current spectrum (1D view)
    const Annotations1DContainer & getCurrentAnnotations() const
    {
      return annotations_1d[current_spectrum_];
    }

    /// Returns a mutable reference to the annotations of the current spectrum (1D view)
    Annotations1DContainer & getCurrentAnnotations()
    {
      return annotations_1d[current_spectrum_];
    }

    /// Returns a const reference to the annotations of the current spectrum (1D view)
    const Annotations1DContainer & getAnnotations(Size spectrum_index) const
    {
      return annotations_1d[spectrum_index];
    }

    /// Returns a mutable reference to the annotations of the current spectrum (1D view)
    Annotations1DContainer & getAnnotations(Size spectrum_index)
    {
      return annotations_1d[spectrum_index];
    }

    /**
    @brief Returns a const reference to the current spectrum (1D view)

    @note Only use this function to access the current spectrum for the 1D view
    */
    const ExperimentType::SpectrumType & getCurrentSpectrum() const;

    void sortCurrentSpectrumByPosition()
    {
      cached_spectrum_.sortByPosition();
    }

    /// Returns a const-copy of the required spectrum which is guaranteed to be populated with raw data
    const ExperimentType::SpectrumType getSpectrum(Size spectrum_idx) const;
      
    /// Get the index of the current spectrum (1D view)
    Size getCurrentSpectrumIndex() const
    {
      return current_spectrum_;
    }

    /// Set the index of the current spectrum (1D view)
    void setCurrentSpectrumIndex(Size index)
    {
      current_spectrum_ = index;
      updateCache_();
    }


    /// get the full chromExperiment
    /// Could be backed up in layer.getChromatogramData() (if layer.getPeakDataMuteable() shows converted chroms already)
    /// ... or layer.getChromatogramData() is empty and thus layer.getPeakDataMuteable() is the original chrom data
    ExperimentSharedPtrType getFullChromData()
    {
      ExperimentSharedPtrType exp_sptr(getChromatogramData().get() == nullptr ||
          getChromatogramData().get()->getNrChromatograms() == 0
             ? getPeakDataMuteable() : getChromatogramData());
      return exp_sptr;
    }

    /// Check whether the current layer should be represented as ion mobility
    bool isIonMobilityData() const
    {
      return this->getPeakData()->size() > 0 &&
             this->getPeakData()->metaValueExists("is_ion_mobility") &&
             this->getPeakData()->getMetaValue("is_ion_mobility").toBool();
    }

    void labelAsIonMobilityData() const
    {
      peak_map_->setMetaValue("is_ion_mobility", "true");
    }

    /// Check whether the current layer contains DIA (SWATH-MS) data
    bool isDIAData() const
    {
      return this->getPeakData()->size() > 0 &&
             this->getPeakData()->metaValueExists("is_dia_data") &&
             this->getPeakData()->getMetaValue("is_dia_data").toBool();
    }

    /// Label the current layer as DIA (SWATH-MS) data
    void labelAsDIAData()
    {
      peak_map_->setMetaValue("is_dia_data", "true");
    }

    /**
    @brief Check whether the current layer is a chromatogram
     
    This is needed because type will *not* distinguish properly between
    chromatogram and spectra data. This is due to the fact that we store 
    chromatograms for display in 1D in a data layer using MSSpectrum and 
    so the layer looks like PEAK data to tools. 
    */
    bool chromatogram_flag_set() const
    {
      return this->getPeakData()->size() > 0 &&
             this->getPeakData()->metaValueExists("is_chromatogram") &&
             this->getPeakData()->getMetaValue("is_chromatogram").toBool();
    }

    /// set the chromatogram flag
    void set_chromatogram_flag()
    {
      peak_map_->setMetaValue("is_chromatogram", "true");
    }

    /// remove the chromatogram flag
    void remove_chromatogram_flag()
    {
      if (this->chromatogram_flag_set())
      {
        peak_map_->removeMetaValue("is_chromatogram");
      }
    }
    
    /**
    @brief Update ranges of all data structures

    Updates ranges of all tracked data structures 
    (spectra, chromatograms, features etc).
    */
    void updateRanges();

    /// Returns the minimum intensity of the internal data, depending on type
    float getMinIntensity() const;

    /// Returns the maximum intensity of the internal data, depending on type
    float getMaxIntensity() const;

    /// updates the PeakAnnotations in the current PeptideHit with manually changed annotations
    /// if no PeptideIdentification or PeptideHit for the spectrum exist, it is generated
    void synchronizePeakAnnotations();

    /// remove peak annotations in the given list from the currently active PeptideHit
    void removePeakAnnotationsFromPeptideHit(const std::vector<Annotation1DItem*>& selected_annotations);

    /// if this layer is visible
    bool visible;

    /// if this layer is flipped (1d mirror view)
    bool flipped;

    /// data type (peak or feature data)
    DataType type;

    private:
    /// layer name
    String name_;

    public:
      const String& getName() const
      {
        return name_;
      }
      void setName(const String& new_name)
      {
        name_ = new_name;
      }

    /// file name of the file the data comes from (if available)
    String filename;

    /// peptide identifications
    std::vector<PeptideIdentification> peptides;

    /// Layer parameters
    Param param;

    /// Gradient for 2D and 3D views
    MultiGradient gradient;

    /// Filters to apply before painting
    DataFilters filters;

    /// Annotations of all spectra of the experiment (1D view)
    std::vector<Annotations1DContainer> annotations_1d;

    /// Peak colors of the currently shown spectrum
    std::vector<QColor> peak_colors_1d;

    /// Flag that indicates if the layer data can be modified (so far used for features only)
    bool modifiable;

    /// Flag that indicates that the layer data was modified since loading it
    bool modified;

    /// Label type
    LabelType label;

    /// Selected peptide id and hit index (-1 if none is selected)
    int peptide_id_index;
    int peptide_hit_index;

    /// get name augmented with attributes, e.g. [flipped], or '*' if modified
    String getDecoratedName() const;

private:
    /// Update current cached spectrum for easy retrieval
    void updateCache_();

    /// updates the PeakAnnotations in the current PeptideHit with manually changed annotations
    void updatePeptideHitAnnotations_(PeptideHit& hit);

    /// feature data
    FeatureMapSharedPtrType features_;

    /// consensus feature data
    ConsensusMapSharedPtrType consensus_map_;

    /// peak data
    ExperimentSharedPtrType peak_map_;

    /// on disc peak data
    ODExperimentSharedPtrType on_disc_peaks;

    /// chromatogram data
    ExperimentSharedPtrType chromatogram_map_;

    /// Chromatogram annotation data
    OSWDataSharedPtrType chrom_annotation_;

    /// Index of the current spectrum
    Size current_spectrum_;

    /// Current cached spectrum
    ExperimentType::SpectrumType cached_spectrum_;
  };

  /// A base class to annotate layers of specific types with (identification) data
  /// @hint Add new derived classes to getAnnotatorWhichSupports() to enable automatic annotation in TOPPView 
  class LayerAnnotatorBase
  {
    public:
      /**
        @brief C'tor with params
        
        @param supported_types Which identification data types are allowed to be opened by the user in annotate()
        @param file_dialog_text The header text of the file dialog shown to the user
        @param gui_lock Optional GUI element which will be locked (disabled) during call to 'annotateWorker_'; can be null_ptr
      **/
      LayerAnnotatorBase(const FileTypes::FileTypeList& supported_types, const String& file_dialog_text, QWidget* gui_lock);
      
      /// Annotates a @p layer, writing messages to @p log and showing QMessageBoxes on errors.
      /// The input file is selected via a file-dialog which is opened with @p current_path as initial path.
      /// The filetype is checked to be one of the supported_types_ before the annotateWorker_ function is called
      /// as implemented by the derived classes
      bool annotateWithFileDialog(LayerData& layer, LogWindow& log, const String& current_path) const;

      /// Annotates a @p layer, given a filename from which to load the data.
      /// The filetype is checked to be one of the supported_types_ before the annotateWorker_ function is called
      /// as implemented by the derived classes
      bool annotateWithFilename(LayerData& layer, LogWindow& log, const String& filename) const;

      /// get a derived annotator class, which supports annotation of the given filetype.
      /// If multiple class support this type (currently not the case) an Exception::IllegalSelfOperation will be thrown
      /// If NO class supports this type, the unique_ptr points to nothing (.get() == nullptr).
      static std::unique_ptr<LayerAnnotatorBase> getAnnotatorWhichSupports(const FileTypes::Type& type);

      /// see getAnnotatorWhichSupports(const FileTypes::Type& type). Filetype is queried from filename
      static std::unique_ptr<LayerAnnotatorBase> getAnnotatorWhichSupports(const String& filename);

    protected:
      /// abstract virtual worker function to annotate a layer using content from the @p filename
      /// returns true on success
      virtual bool annotateWorker_(LayerData& layer, const String& filename, LogWindow& log) const = 0;
      
      const FileTypes::FileTypeList supported_types_;
      const String file_dialog_text_;
      QWidget* gui_lock_ = nullptr; ///< optional widget which will be locked when calling annotateWorker_() in child-classes
  };

  /// Annotate a layer with PeptideIdentifications using Layer::annotate(pepIDs, protIDs).
  /// The ID data is loaded from a file selected by the user via a file-dialog.
  class LayerAnnotatorPeptideID
    : public LayerAnnotatorBase
  {
    public:
      LayerAnnotatorPeptideID(QWidget* gui_lock)
       : LayerAnnotatorBase(std::vector<FileTypes::Type>{ FileTypes::IDXML, FileTypes::MZIDENTML },
                            "Select peptide identification data", gui_lock)
      {}

  protected:
    /// loads the ID data from @p filename and calls Layer::annotate.
    /// Always returns true (unless an exception is thrown from internal sub-functions)
    virtual bool annotateWorker_(LayerData& layer, const String& filename, LogWindow& log) const;
  };

  /// Annotate a layer with AccurateMassSearch results (from an AMS-featureXML file).
  /// The featuremap is loaded from a file selected by the user via a file-dialog.
  class LayerAnnotatorAMS
    : public LayerAnnotatorBase
  {
  public:
    LayerAnnotatorAMS(QWidget* gui_lock)
      : LayerAnnotatorBase(std::vector<FileTypes::Type>{ FileTypes::FEATUREXML },
                           "Select AccurateMassSearch's featureXML file", gui_lock)
    {}

  protected:
    /// loads the featuremap from @p filename and calls Layer::annotate.
    /// Returns false if featureXML file was not created by AMS, and true otherwise (unless an exception is thrown from internal sub-functions)
    virtual bool annotateWorker_(LayerData& layer, const String& filename, LogWindow& log) const;
  };
  
  /// Annotate a chromatogram layer with ID data (from an OSW sqlite file as produced by OpenSwathWorkflow or pyProphet).
  /// The OSWData is loaded from a file selected by the user via a file-dialog.
  class LayerAnnotatorOSW
    : public LayerAnnotatorBase
  {
  public:
    LayerAnnotatorOSW(QWidget* gui_lock)
      : LayerAnnotatorBase(std::vector<FileTypes::Type>{ FileTypes::OSW },
                           "Select OpenSwath/pyProphet output file", gui_lock)
    {}

  protected:
    /// loads the OSWData from @p filename and stores the data using Layer::setChromatogramAnnotation()
    /// Always returns true (unless an exception is thrown from internal sub-functions)
    virtual bool annotateWorker_(LayerData& layer, const String& filename, LogWindow& log) const;
  };

  /// Print the contents to a stream.
  OPENMS_GUI_DLLAPI std::ostream& operator<<(std::ostream & os, const LayerData & rhs);

} //namespace


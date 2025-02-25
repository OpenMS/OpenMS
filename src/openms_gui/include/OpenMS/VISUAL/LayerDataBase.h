// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/PROCESSING/MISC/DataFilters.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/VISUAL/LogWindow.h>
#include <OpenMS/VISUAL/MISC/CommonDefs.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <boost/shared_ptr.hpp>

#include <bitset>
#include <vector>

class QWidget;

namespace OpenMS
{
  class LayerData1DBase;
  class LayerStoreData;
  class LayerStatistics;
  class OSWData;
  class Painter1DBase;
  class Painter2DBase;

  template <int N_DIM> class DimMapper;

  struct LayerDataDefs
  {

    /// Result of computing a projection on X and Y axis in a 2D Canvas; see LayerDataBase::getProjection()
    struct ProjectionData
    {
      /// C'tor
      ProjectionData();
      /// Move C'tor
      ProjectionData(ProjectionData&&);
      /// D'tor
      ~ProjectionData(); // needs to be implemented in cpp, since inline would require a full definition of LayerData1DBase;

      std::unique_ptr<LayerData1DBase> projection_ontoX;
      std::unique_ptr<LayerData1DBase> projection_ontoY;
      struct Summary {
        UInt number_of_datapoints {0};
        Peak1D::IntensityType max_intensity {0};
        double sum_intensity {0}; // double since sum could get large
      } stats;
    };

    /** @name Type definitions */
    //@{
    /// Dataset types.
    /// Order in the enum determines the order in which layer types are drawn.
    enum DataType
    {
      DT_PEAK,         ///< Spectrum profile or centroided data
      DT_CHROMATOGRAM, ///< Chromatogram data
      DT_FEATURE,      ///< Feature data
      DT_CONSENSUS,    ///< Consensus feature data
      DT_IDENT,        ///< Peptide identification data
      DT_UNKNOWN       ///< Undefined data type indicating an error
    };

    /// Flags that determine which information is shown.
    enum Flags
    {
      F_HULL,        ///< Features: Overall convex hull
      F_HULLS,       ///< Features: Convex hulls of single mass traces
      F_UNASSIGNED,  ///< Features: Unassigned peptide hits
      P_PRECURSORS,  ///< Peaks: Mark precursor peaks of MS/MS scans
      P_PROJECTIONS, ///< Peaks: Show projections
      C_ELEMENTS,    ///< Consensus features: Show elements
      I_PEPTIDEMZ,   ///< Identifications: m/z source
      I_LABELS,      ///< Identifications: Show labels (not sequences)
      SIZE_OF_FLAGS
    };

    /// Label used in visualization
    enum LabelType
    {
      L_NONE,       ///< No label is displayed
      L_INDEX,      ///< The element number is used
      L_META_LABEL, ///< The 'label' meta information is used
      L_ID,         ///< The best peptide hit of the first identification run is used
      L_ID_ALL,     ///< All peptide hits of the first identification run are used
      SIZE_OF_LABEL_TYPE
    };

    /// Label names
    static const std::string NamesOfLabelType[SIZE_OF_LABEL_TYPE];

    /// Features
    typedef FeatureMap FeatureMapType;

    /// SharedPtr on feature map
    typedef boost::shared_ptr<FeatureMap> FeatureMapSharedPtrType;

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
  };

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
  a vector of LayerDataBase objects.

  @ingroup PlotWidgets
  */
#ifdef _MSC_VER
  #pragma warning(disable : 4250) // 'class1' : inherits 'class2::member' via dominance
#endif
  class OPENMS_GUI_DLLAPI LayerDataBase : public LayerDataDefs
  {
  public:
    /// Actual state of each flag
    std::bitset<SIZE_OF_FLAGS> flags;

    //@}

    /// Default constructor (for virtual inheritance)
    LayerDataBase() = delete;     // <-- this is the problem. Call assignment op in 1DPeak???
    /// C'tor for child classes
    explicit LayerDataBase(const DataType type) : type(type) {}
    /// Copy-C'tor
    LayerDataBase(const LayerDataBase& ld) = default;
    /// Assignment operator
    LayerDataBase& operator=(const LayerDataBase& ld) = delete;
    /// Move-C'tor - do not move from this class since its a virtual base class (diamond problem) and the move c'tor may be called twice (which would loose data!)
    /// Instead of painstakingly writing user-defined move c'tors which check for moving for all the direct child classes, 
    /// we'd rather use copy (which is the automatic fallback, and safe) and incur a small performance hit
    LayerDataBase(LayerDataBase&& ld) = delete;
    /// Move assignment -- deleted, by same argument as for move c'tor
    LayerDataBase& operator=(LayerDataBase&& ld) = delete;
    /// D'tor
    virtual ~LayerDataBase() = default;

    /**
     * \brief Obtain a painter which can draw the layer on a 2D canvas
     * \return A painter
     */
    virtual std::unique_ptr<Painter2DBase> getPainter2D() const = 0;


    /**
     * \brief Create a shallow copy (i.e. shared experimental data using shared_ptr) of the current layer, and make it 1D (i.e. support showing a single spec/chrom etc)
     * \return A new layer for 1D
     */
    virtual std::unique_ptr <LayerData1DBase> to1DLayer() const = 0;

    /// Returns a visitor which contains the current visible data and can write the data to disk
    virtual std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& /*visible_range*/, const DataFilters& /*layer_filters*/) const
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    /// Returns a visitor which contains the the full data of the layer and can write the data to disk in the appropriate format (e.g. mzML)
    virtual std::unique_ptr<LayerStoreData> storeFullData() const
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /// Calculate a projection of the current layer for the given unit and the given area.
    /// E.g. the area might be restricted in RT and m/z, and then requested projection should return the XIC (@p unit == RT)
    /// It is up to the implementation to decide on binning.
    // todo: put this into a LayerData2DBase class, since a LayerData1DPeak should not implement this.
    virtual ProjectionData getProjection(const DIM_UNIT unit_x, const DIM_UNIT unit_y, const RangeAllType& area) const = 0;

    /**
     * \brief Find the closest datapoint within the given range and return a proxy to that datapoint
     * \param area Range to search in. Only dimensions used in the canvas are populated.
     * \return A proxy (e.g. scan + peak index in an MSExperiment) which points to the data
     */
    virtual PeakIndex findClosestDataPoint(const RangeAllType& area) const
    {
      (void)area; // allow doxygen to document the param
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /**
     * \brief Find the datapoint with the highest intensity within the given range and return a proxy to that datapoint
     * \param area Range to search in. Only dimensions used in the canvas are populated.
     * \return A proxy (e.g. scan + peak index in an MSExperiment) which points to the data
     */
    virtual PeakIndex findHighestDataPoint(const RangeAllType& area) const
    {
      (void)area; // allow doxygen to document the param
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }


    /**
     * \brief Convert a PeakIndex to a XY coordinate (via @p mapper).
     * \param peak The Peak to convert
     * \param mapper Converts the internal representation (e.g. Peak1D) to an XY coordinate
     * \return XY coordinate in data units (e.g. X=m/z, Y=intensity)
     */
    virtual PointXYType peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const = 0;

    /**
     * \brief Get name and value of all data-arrays corresponding to the given datapoint
     *
     * Empty (or shorter) data-arrays are skipped.
     *
     * \param peak_index The datapoint
     * \return A string, e.g. "fwhm: 20, im: 3.3", depending on which float/string dataarrays are populated for the given datapoint
     */
    virtual String getDataArrayDescription(const PeakIndex& peak_index)
    {
      (void)peak_index; // allow doxygen to document the param
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /// add peptide identifications to the layer
    /// Only supported for DT_PEAK, DT_FEATURE and DT_CONSENSUS.
    /// Will return false otherwise.
    virtual bool annotate(const std::vector<PeptideIdentification>& identifications,
                          const std::vector<ProteinIdentification>& protein_identifications);


    /**
      @brief Update ranges of the underlying data
    */
    virtual void updateRanges() = 0;

    /// Returns the minimum intensity of the internal data, depending on type
    float getMinIntensity() const;

    /// Returns the maximum intensity of the internal data, depending on type
    float getMaxIntensity() const;

    using RangeAllType = RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>;

    /// Returns the data range of the whole layer (i.e. all scans/chroms/etc) in all known dimensions.
    /// If a layer does not support the dimension (or the layer is empty) the dimension will be empty
    /// If you need the data range for a 1D view (i.e. only a single spec/chrom/etc), call 'LayerDataBase1D::getRange1D()'
    virtual RangeAllType getRange() const = 0;

    /// Compute layer statistics (via visitor)
    virtual std::unique_ptr<LayerStatistics> getStats() const = 0;

    /// The name of the layer, usually the basename of the file
    const String& getName() const
    {
      return name_;
    }
    /// Set the name of the layer, usually the basename of the file
    void setName(const String& new_name)
    {
      name_ = new_name;
    }

    /// get the extra annotation to the layers name, e.g. '[39]' for which chromatogram index is currently shown in 1D
    const String& getNameSuffix() const
    {
      return name_suffix_;
    }
    /// set an extra annotation as suffix to the layers name, e.g. '[39]' for which chromatogram index is currently shown in 1D
    void setNameSuffix(const String& decorator)
    {
      name_suffix_ = decorator;
    }

    /// get name augmented with attributes, e.g. '*' if modified
    virtual String getDecoratedName() const;

    /// if this layer is visible
    bool visible = true;

    /// data type (peak or feature data, etc)
    DataType type = DT_UNKNOWN;

    /// file name of the file the data comes from (if available)
    String filename;

    /// Layer parameters
    Param param;

    /// Gradient for 2D and 3D views
    MultiGradient gradient;

    /// Filters to apply before painting
    DataFilters filters;

    /// Flag that indicates if the layer data can be modified (so far used for features only)
    bool modifiable = false;

    /// Flag that indicates that the layer data was modified since loading it
    bool modified = false;

    /// Label type
    LabelType label = L_NONE;

    /// Selected peptide id and hit index (-1 if none is selected)
    int peptide_id_index = -1;
    int peptide_hit_index = -1;

  private:
    /// layer name
    String name_;
    /// an extra annotation as suffix to the layers name, e.g. '[39]' for which chromatogram index is currently shown in 1D
    String name_suffix_;
  };

  /// A base class to annotate layers of specific types with (identification) data
  /// 
  /// @note Add new derived classes to getAnnotatorWhichSupports() to enable automatic annotation in TOPPView
  class LayerAnnotatorBase
  {
  public:
    /**
        @brief C'tor with params
        
        @param supported_types Which identification data types are allowed to be opened by the user in annotate()
        @param file_dialog_text The header text of the file dialog shown to the user
        @param gui_lock Optional GUI element which will be locked (disabled) during call to 'annotateWorker_'; can be null_ptr
      **/
    LayerAnnotatorBase(const FileTypeList& supported_types, const String& file_dialog_text, QWidget* gui_lock);
    
    /// Make D'tor virtual for correct destruction from pointers to base
    virtual ~LayerAnnotatorBase() = default;

    /// Annotates a @p layer, writing messages to @p log and showing QMessageBoxes on errors.
    /// The input file is selected via a file-dialog which is opened with @p current_path as initial path.
    /// The file type is checked to be one of the supported_types_ before the annotateWorker_ function is called
    /// as implemented by the derived classes
    bool annotateWithFileDialog(LayerDataBase& layer, LogWindow& log, const String& current_path) const;

    /// Annotates a @p layer, given a filename from which to load the data.
    /// The file type is checked to be one of the supported_types_ before the annotateWorker_ function is called
    /// as implemented by the derived classes
    bool annotateWithFilename(LayerDataBase& layer, LogWindow& log, const String& filename) const;

    /// get a derived annotator class, which supports annotation of the given file type.
    /// If multiple class support this type (currently not the case) an Exception::IllegalSelfOperation will be thrown
    /// If NO class supports this type, the unique_ptr points to nothing (.get() == nullptr).
    static std::unique_ptr<LayerAnnotatorBase> getAnnotatorWhichSupports(const FileTypes::Type& type);

    /// see getAnnotatorWhichSupports(const FileTypes::Type& type). File type is queried from filename
    static std::unique_ptr<LayerAnnotatorBase> getAnnotatorWhichSupports(const String& filename);

  protected:
    /// abstract virtual worker function to annotate a layer using content from the @p filename
    /// returns true on success
    virtual bool annotateWorker_(LayerDataBase& layer, const String& filename, LogWindow& log) const = 0;

    const FileTypeList supported_types_;
    const String file_dialog_text_;
    QWidget* gui_lock_ = nullptr;///< optional widget which will be locked when calling annotateWorker_() in child-classes
  };

  /// Annotate a layer with PeptideIdentifications using Layer::annotate(pepIDs, protIDs).
  /// The ID data is loaded from a file selected by the user via a file-dialog.
  class LayerAnnotatorPeptideID : public LayerAnnotatorBase
  {
  public:
    LayerAnnotatorPeptideID(QWidget* gui_lock) :
        LayerAnnotatorBase(std::vector<FileTypes::Type>{FileTypes::IDXML, FileTypes::MZIDENTML},
                           "Select peptide identification data", gui_lock)
    {
    }

  protected:
    /// loads the ID data from @p filename and calls Layer::annotate.
    /// Always returns true (unless an exception is thrown from internal sub-functions)
    virtual bool annotateWorker_(LayerDataBase& layer, const String& filename, LogWindow& log) const;
  };

  /// Annotate a layer with AccurateMassSearch results (from an AMS-featureXML file).
  /// The featuremap is loaded from a file selected by the user via a file-dialog.
  class LayerAnnotatorAMS : public LayerAnnotatorBase
  {
  public:
    LayerAnnotatorAMS(QWidget* gui_lock) :
        LayerAnnotatorBase(std::vector<FileTypes::Type>{FileTypes::FEATUREXML},
                           "Select AccurateMassSearch's featureXML file", gui_lock)
    {
    }

  protected:
    /// loads the featuremap from @p filename and calls Layer::annotate.
    /// Returns false if featureXML file was not created by AMS, and true otherwise (unless an exception is thrown from internal sub-functions)
    virtual bool annotateWorker_(LayerDataBase& layer, const String& filename, LogWindow& log) const;
  };

  /// Annotate a chromatogram layer with ID data (from an OSW sqlite file as produced by OpenSwathWorkflow or pyProphet).
  /// The OSWData is loaded from a file selected by the user via a file-dialog.
  class LayerAnnotatorOSW : public LayerAnnotatorBase
  {
  public:
    LayerAnnotatorOSW(QWidget* gui_lock) :
        LayerAnnotatorBase(std::vector<FileTypes::Type>{FileTypes::OSW},
                           "Select OpenSwath/pyProphet output file", gui_lock)
    {
    }

  protected:
    /// loads the OSWData from @p filename and stores the data using Layer::setChromatogramAnnotation()
    /// Always returns true (unless an exception is thrown from internal sub-functions)
    virtual bool annotateWorker_(LayerDataBase& layer, const String& filename, LogWindow& log) const;
  };

  /// Print the contents to a stream.
  OPENMS_GUI_DLLAPI std::ostream& operator<<(std::ostream& os, const LayerDataBase& rhs);

}// namespace OpenMS


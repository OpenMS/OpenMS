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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>
#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>

#include <boost/shared_ptr.hpp>

#include <vector>
#include <bitset>

namespace OpenMS
{
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

  @note Layer is mainly used as a member variable of SpectrumCanvas which holds
  a vector of LayerData objects.

  @ingroup SpectrumWidgets
  */
  class LayerData
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

    //@}

    /// Default constructor
    LayerData() :
      flags(),
      visible(true),
      flipped(false),
      type(DT_UNKNOWN),
      name(),
      filename(),
      peptides(),
      param(),
      gradient(),
      filters(),
      annotations_1d(),
      peak_colors_1d(),
      modifiable(false),
      modified(false),
      label(L_NONE),
      peptide_id_index(-1),
      peptide_hit_index(-1),
      features(new FeatureMapType()),
      consensus(new ConsensusMapType()),
      peaks(new ExperimentType()),
      on_disc_peaks(new OnDiscMSExperiment()),
      chromatograms(new ExperimentType()),
      current_spectrum_(0),
      cached_spectrum_()
    {
      annotations_1d.resize(1);
    }

    /// Returns a const reference to the current feature data
    const FeatureMapSharedPtrType & getFeatureMap() const
    {
      return features;
    }

    /// Returns a const reference to the current feature data
    FeatureMapSharedPtrType & getFeatureMap()
    {
      return features;
    }

    /// Returns a const reference to the consensus feature data
    const ConsensusMapSharedPtrType & getConsensusMap() const
    {
      return consensus;
    }

    /// Returns current consensus map (mutable)
    ConsensusMapSharedPtrType & getConsensusMap()
    {
      return consensus;
    }

    /**
    @brief Returns a const reference to the current in-memory peak data

    @note Depending on the caching strategy (on-disk or in-memory), all or some
    spectra may have zero size and contain only meta data since peak data is
    cached on disk.

    @note Do *not* use this function to access the current spectrum for the 1D view
    */
    const ConstExperimentSharedPtrType getPeakData() const;

    /**
    @brief Returns a mutable reference to the current in-memory peak data

    @note Depending on the caching strategy (on-disk or in-memory), all or some
    spectra may have zero size and contain only meta data since peak data is
    cached on disk.

    @note Do *not* use this function to access the current spectrum for the 1D view
    */
    const ExperimentSharedPtrType & getPeakDataMuteable() {return peaks;}

    /**
    @brief Set the current in-memory peak data
    */
    void setPeakData(ExperimentSharedPtrType p)
    {
      peaks = p;
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
      return chromatograms;
    }

    /// Returns a mutable reference to the current chromatogram data
    ExperimentSharedPtrType & getChromatogramData()
    {
      return chromatograms;
    }

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
    const ExperimentType::SpectrumType getSpectrum(Size spectrum_idx) const
    {
      if (spectrum_idx == current_spectrum_) return cached_spectrum_;

      if ((*peaks)[spectrum_idx].size() > 0)
      {
        return (*peaks)[spectrum_idx];
      }
      else if (!on_disc_peaks->empty())
      {
        return on_disc_peaks->getSpectrum(spectrum_idx);
      }
      return (*peaks)[spectrum_idx];
    }
      
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

    /// Check whether the current layer should be represented as ion mobility
    bool isIonMobilityData() const
    {
      return this->getPeakData()->size() > 0 &&
             this->getPeakData()->metaValueExists("is_ion_mobility") &&
             this->getPeakData()->getMetaValue("is_ion_mobility").toBool();
    }

    void labelAsIonMobilityData() const
    {
      peaks->setMetaValue("is_ion_mobility", "true");
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
      peaks->setMetaValue("is_dia_data", "true");
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
      peaks->setMetaValue("is_chromatogram", "true");
    }

    /// remove the chromatogram flag
    void remove_chromatogram_flag()
    {
      if (this->chromatogram_flag_set())
      {
        peaks->removeMetaValue("is_chromatogram");
      }
    }
    
    /**
    @brief Update ranges of all data structures

    Updates ranges of all tracked data structures 
    (spectra, chromatograms, features etc).
    */
    void updateRanges();

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

    /// layer name
    String name;

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
    String getDecoratedName() const
    {
      String n = name;
      if (flipped)
      {
        n += " [flipped]";
      }
      if (modified)
      {
        n += '*';
      }
      return n;
    }

private:


    /// Update current cached spectrum for easy retrieval
    void updateCache_();

    /// updates the PeakAnnotations in the current PeptideHit with manually changed annotations
    void updatePeptideHitAnnotations_(PeptideHit& hit);

    /// feature data
    FeatureMapSharedPtrType features;

    /// consensus feature data
    ConsensusMapSharedPtrType consensus;

    /// peak data
    ExperimentSharedPtrType peaks;

    /// on disc peak data
    ODExperimentSharedPtrType on_disc_peaks;

    /// chromatogram data
    ExperimentSharedPtrType chromatograms;

    /// Index of the current spectrum
    Size current_spectrum_;

    /// Current cached spectrum
    ExperimentType::SpectrumType cached_spectrum_;

  };

  /// Print the contents to a stream.
  OPENMS_GUI_DLLAPI std::ostream & operator<<(std::ostream & os, const LayerData & rhs);

} //namespace


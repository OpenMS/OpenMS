// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_VISUAL_LAYERDATA_H
#define OPENMS_VISUAL_LAYERDATA_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
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

      @ingroup SpectrumWidgets
  */
  class LayerData
  {
public:
    /** @name Type definitions */
    //@{
    /// Dataset types
    enum DataType
    {
      DT_PEAK,                ///< Spectrum profile or centroided data
      DT_FEATURE,         ///< Feature data
      DT_CONSENSUS,       ///< Consensus feature data
      DT_CHROMATOGRAM,    ///< Chromatogram data
      DT_IDENT,           ///< Peptide identification data
      DT_UNKNOWN              ///< Undefined data type indicating an error
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
      SIZE_OF_FLAGS
    };

    /// Actual state of each flag
    std::bitset<SIZE_OF_FLAGS> flags;

    /// Label used in visualization
    enum LabelType
    {
      L_NONE,                           ///< No label is displayed
      L_INDEX,                          ///< The element number is used
      L_META_LABEL,                 ///< The 'label' meta information is used
      L_ID,                                 ///< The best peptide hit of the first identification run is used
      L_ID_ALL,                         ///< All peptide hits of the first identification run are used
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
      modifiable(false),
      modified(false),
      label(L_NONE),
      peptide_id_index(-1),
      peptide_hit_index(-1),
      features(new FeatureMapType()),
      consensus(new ConsensusMapType()),
      peaks(new ExperimentType()),
      chromatograms(new ExperimentType()),
      current_spectrum_(0)
    {
      annotations_1d.resize(1);
    }

    /// Returns a const reference to the current spectrum (1d view)
    const ExperimentType::SpectrumType & getCurrentSpectrum() const;

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

    /// Returns a const reference to the current peak data
    const ExperimentSharedPtrType & getPeakData() const
    {
      return peaks;
    }

    /// Returns a mutable reference to the current peak data
    ExperimentSharedPtrType & getPeakData()
    {
      return peaks;
    }

    /// Returns a const reference to the current chromatogram data
    const ExperimentSharedPtrType & getChromatogramData() const
    {
      return chromatograms;
    }

    /// Returns a mutable reference to the current chromatogram data
    ExperimentSharedPtrType & getChromatogramData()
    {
      return chromatograms;
    }

    /// Returns a const reference to the annotations of the current spectrum (1d view)
    const Annotations1DContainer & getCurrentAnnotations() const
    {
      return annotations_1d[current_spectrum_];
    }

    /// Returns a mutable reference to the annotations of the current spectrum (1d view)
    Annotations1DContainer & getCurrentAnnotations()
    {
      return annotations_1d[current_spectrum_];
    }

    /// Returns a const reference to the annotations of the current spectrum (1d view)
    const Annotations1DContainer & getAnnotations(Size spectrum_index) const
    {
      return annotations_1d[spectrum_index];
    }

    /// Returns a mutable reference to the annotations of the current spectrum (1d view)
    Annotations1DContainer & getAnnotations(Size spectrum_index)
    {
      return annotations_1d[spectrum_index];
    }

    /// Returns a mutable reference to the current spectrum (1d view)
    ExperimentType::SpectrumType & getCurrentSpectrum()
    {
      return (*peaks)[current_spectrum_];
    }

    /// Get the index of the current spectrum
    Size getCurrentSpectrumIndex() const
    {
      return current_spectrum_;
    }

    /// Set the index of the current spectrum
    void setCurrentSpectrumIndex(Size index)
    {
      current_spectrum_ = index;
    }

    /// Check whether the current layer is a chromatogram
    // we need this specifically because this->type will *not* distinguish
    // chromatogram and spectra data since we need to store chromatograms for
    // the 1D case in a layer that looks like PEAK data to all tools.
    bool chromatogram_flag_set() const
    {
      return this->getPeakData()->size() > 0 &&
             this->getPeakData()->metaValueExists("is_chromatogram") &&
             this->getPeakData()->getMetaValue("is_chromatogram").toBool();
    }

    // set the chromatogram flag
    void set_chromatogram_flag()
    {
      this->getPeakData()->setMetaValue("is_chromatogram", "true");
    }

    // remove the chromatogram flag
    void remove_chromatogram_flag()
    {
      if (this->chromatogram_flag_set())
      {
        this->getPeakData()->removeMetaValue("is_chromatogram");
      }
    }

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

    /// Flag that indicates if the layer data can be modified (so far used for features only)
    bool modifiable;

    /// Flag that indicates that the layer data was modified since loading it
    bool modified;

    /// Label type
    LabelType label;

    /// Selected peptide id and hit index (-1 if none is selected)
    int peptide_id_index;
    int peptide_hit_index;

private:
    /// updates the PeakAnnotations in the current PeptideHit with manually changed annotations
    void updatePeptideHitAnnotations_(PeptideHit& hit);

    /// feature data
    FeatureMapSharedPtrType features;

    /// consensus feature data
    ConsensusMapSharedPtrType consensus;

    /// peak data
    ExperimentSharedPtrType peaks;

    /// chromatogram data
    ExperimentSharedPtrType chromatograms;

    /// Index of the current spectrum
    Size current_spectrum_;
  };

  /// Print the contents to a stream.
  OPENMS_GUI_DLLAPI std::ostream & operator<<(std::ostream & os, const LayerData & rhs);

} //namespace

#endif

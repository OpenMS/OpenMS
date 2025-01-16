// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/FeatureMap.h>


#include <vector>
#include <list>
#include <map>

namespace OpenMS
{
  typedef Peak2D PeakType;
  class string;

  /** @brief A container type that gathers peaks similar in m/z and moving along retention time.

    Depending on the method of extraction a mass trace could virtually
    represent a complete extracted ion chromatogram (XIC) or merely a part of
    it (e.g., a chromatographic peak). The kernel class provides methods for
    computing mass trace characteristics such as its centroid m/z and retention
    time. Coeluting mass traces can be further assembled to complete isotope
    patterns of peptides/metabolites.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI MassTrace
  {
public:

    // must match to names_of_quantmethod[]
    enum MT_QUANTMETHOD {
      MT_QUANT_AREA = 0,  ///< quantify by area
      MT_QUANT_MEDIAN,    ///< quantify by median of intensities
      MT_QUANT_HEIGHT,    ///< quantify by peak height
      SIZE_OF_MT_QUANTMETHOD
    };
    static const std::string names_of_quantmethod[SIZE_OF_MT_QUANTMETHOD];

    /// converts a string to enum value; returns 'SIZE_OF_MT_QUANTMETHOD' upon error
    static MT_QUANTMETHOD getQuantMethod(const String& val);

    /** @name Constructors and Destructor
    */
    ///@{

    /// Default constructor
    MassTrace() = default;

    /// Detailed constructor 
    /// (useful, since Mass Traces are commonly assembled by prepending and appending -- which is faster using lists)
    MassTrace(const std::list<PeakType>& trace_peaks);

    /// Detailed constructor for vector
    MassTrace(const std::vector<PeakType>& trace_peaks);

    /// Destructor
    ~MassTrace() = default;

    /// Copy constructor
    MassTrace(const MassTrace &) = default;

    /// Assignment operator
    MassTrace & operator=(const MassTrace &) = default;

    /// Random access operator
    PeakType& operator[](const Size & mt_idx);
    const PeakType& operator[](const Size & mt_idx) const;
    ///@}

    /** @name Iterators
        @brief Enables mutable/immutable access to the mass trace's peaks.
    */
    ///@{
    typedef std::vector<PeakType>::iterator iterator;
    typedef std::vector<PeakType>::const_iterator const_iterator;
    typedef std::vector<PeakType>::reverse_iterator reverse_iterator;
    typedef std::vector<PeakType>::const_reverse_iterator const_reverse_iterator;

    iterator begin()
    {
      return trace_peaks_.begin();
    }

    iterator end()
    {
      return trace_peaks_.end();
    }

    const_iterator begin() const
    {
      return trace_peaks_.begin();
    }

    const_iterator end() const
    {
      return trace_peaks_.end();
    }

    reverse_iterator rbegin()
    {
      return trace_peaks_.rbegin();
    }

    reverse_iterator rend()
    {
      return trace_peaks_.rend();
    }

    const_reverse_iterator rbegin() const
    {
      return trace_peaks_.rbegin();
    }

    const_reverse_iterator rend() const
    {
      return trace_peaks_.rend();
    }
    ///@}

    /** @name Accessor methods
    */

    ///@{

    /// Returns the number of peaks contained in the mass trace.
    Size getSize() const
    {
      return trace_peaks_.size();
    }

    /// Gets label of mass trace.
    String getLabel() const
    {
      return label_;
    }

    /// Sets label of mass trace.
    void setLabel(const String & label)
    {
      label_ = label;
    }

    /// Returns the centroid m/z.
    double getCentroidMZ() const
    {
      return centroid_mz_;
    }

    /// Returns the centroid RT.
    double getCentroidRT() const
    {
      return centroid_rt_;
    }

    double getCentroidSD() const
    {
      return centroid_sd_;
    }

    void setCentroidSD(const double & tmp_sd)
    {
      centroid_sd_ = tmp_sd;
    }

    double getFWHM() const
    {
        return fwhm_;
    }

    /// Returns the length of the trace (as difference in RT)
    double getTraceLength() const
    {
        double length(0.0);

        if (trace_peaks_.size() > 1)
        {
            length = std::fabs(trace_peaks_.rbegin()->getRT() - trace_peaks_.begin()->getRT());
        }

        return length;
    }

    std::pair<Size, Size> getFWHMborders() const
    {
      return std::make_pair(fwhm_start_idx_, fwhm_end_idx_);
    }

    /// Gets smoothed intensities (empty if no smoothing was explicitly done beforehand!).
    const std::vector<double>& getSmoothedIntensities() const
    {
      return smoothed_intensities_;
    }

    /// Set smoothed intensities (smoothing is done externally, e.g. by LowessSmoothing).
    void setSmoothedIntensities(const std::vector<double> & db_vec)
    {
      if (trace_peaks_.size() != db_vec.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Number of smoothed intensities deviates from mass trace size! Aborting...", String(db_vec.size()));
      }

      smoothed_intensities_ = db_vec;
    }

    /// Get average scan time of mass trace
    double getAverageMS1CycleTime() const
    {
      if (trace_peaks_.size() <= 1) return 0.0;

      return (trace_peaks_.rbegin()->getRT() - trace_peaks_.begin()->getRT()) / (trace_peaks_.size() - 1);
    }
    ///@}

    /** @name Computational methods
    */
    ///@{

    /// Sum all non-negative (smoothed!) intensities in the mass trace
    double computeSmoothedPeakArea() const;

    /// Compute area of peaks in the mass trace
    double computePeakArea() const;

    /// Sum all peak intensities in the mass trace
    double computeIntensitySum() const;

    /// Return the index of the mass trace's highest peak within the MassTrace container (based either on raw or smoothed intensities).
    Size findMaxByIntPeak(bool use_smoothed_ints = false) const;

    /// Estimate FWHM of chromatographic peak in seconds (based on either raw or smoothed intensities).
    /// Uses linear interpolation of the two closest points to the half_max intensity in order to get the RT values at exactly the half_max
    /// stores result internally, use getFWHM().
    double estimateFWHM(bool use_smoothed_ints = false);

    /// determine if area or median is used for quantification
    void setQuantMethod(MT_QUANTMETHOD method);

    /// check if area or median is used for quantification
    MT_QUANTMETHOD getQuantMethod() const;

    /// Compute chromatographic peak area within the FWHM range.
    double computeFwhmAreaSmooth() const;
    double computeFwhmArea() const;
    // double computeFwhmAreaSmoothRobust() const;
    // double computeFwhmAreaRobust() const;

    double getIntensity(bool smoothed) const;
    double getMaxIntensity(bool smoothed) const;

    /// Return the mass trace's convex hull.
    ConvexHull2D getConvexhull() const;
    ///@}

    /** @name Update methods for centroid RT and m/z
    */
    ///@{

    void updateSmoothedMaxRT();

    /// Compute & update centroid RT as a intensity-weighted mean of RTs.
    void updateWeightedMeanRT();

    void updateSmoothedWeightedMeanRT();

    /// Compute & update centroid RT as median position of intensities.
    void updateMedianRT();

    /// Compute & update centroid m/z as median of m/z values.
    void updateMedianMZ();

    /// Compute & update centroid m/z as mean of m/z values.
    void updateMeanMZ();

    /// Compute & update centroid m/z as weighted mean of m/z values.
    void updateWeightedMeanMZ();

    /** @brief Compute & update m/z standard deviation of mass trace as weighted mean of m/z values.
    
      Make sure to call update(Weighted)(Mean|Median)MZ() first! <br>
      use getCentroidSD() to get result
    */
    void updateWeightedMZsd();
    ///@}

    /// Average FWHM of m/z peaks
    double fwhm_mz_avg = 0;

private:

    /// median of trace intensities
    double computeMedianIntensity_() const;

    /// calculate x coordinate of start/end indexes at half_max
    /// calculation is based on (yB - yA) / (xB - xA) = (y_eval - yA) / (xC - xA)
    /// solve for xC: xC = xA + ((y_eval - yA) * (xB - xA) / (yB - yA))
    double linearInterpolationAtY_(double xA, double xB, double yA, double yB, double y_eval) const;

    /// Actual MassTrace container for doing centroid calculation, peak width estimation etc.
    std::vector<PeakType> trace_peaks_;

    /// Centroid m/z
    double centroid_mz_ = 0.0;

    /// intensity-weighted STD
    double centroid_sd_ = 0.0;

    /// Centroid RT
    double centroid_rt_ = 0.0;

    /// Trace label
    String label_;

    /// Container for smoothed intensities. Smoothing must be done externally.
    std::vector<double> smoothed_intensities_;

    double fwhm_ = 0.0; ///< FWHM of RT peak
    Size fwhm_start_idx_ = 0; ///< index into 'trace_peaks_' vector (inclusive)
    Size fwhm_end_idx_ = 0; ///< index into 'trace_peaks_' vector (inclusive)

    /// use area under mass trace or the median of intensities
    MT_QUANTMETHOD quant_method_ = MT_QUANT_AREA;
    
  };

}


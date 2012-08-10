// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MASSTRACE_H
#define OPENMS_KERNEL_MASSTRACE_H

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/FeatureMap.h>


#include <vector>
#include <list>
#include <map>


namespace OpenMS
{
typedef Peak2D PeakType;

/** @brief A container type that gathers peaks similar in m/z and moving along retention time.

        Depending on the method of extraction a mass trace could virtually represent a complete ion chromatogram (XIC) or merely a part of it
    (e.g., a chromatographic peak). The kernel class provides methods for computing mass trace characteristics such
      as its centroid m/z and retention time. Coeluting mass traces can be further assembled to complete isotope patterns of peptides/metabolites.

                @ingroup Kernel
        */
class OPENMS_DLLAPI MassTrace
{
public:
    /** @name Constructors and Destructor
        */
    /// Default constructor
    MassTrace();

    /// Detailed constructor 1
    MassTrace(const std::list<PeakType>&, const DoubleReal& scan_time = 1.0);

    /// Detailed constructor 2
    MassTrace(const std::vector<PeakType>&, const DoubleReal& scan_time = 1.0);

    /// Destructor
    ~MassTrace();

    /// Copy constructor
    MassTrace(const MassTrace&);

    /// Assignment operator
    MassTrace& operator= (const MassTrace&);

    /** @name Iterators
      @brief Enables mutable/immutable access to the mass trace's peaks.
        */
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

    /** @name Accessor methods
        */

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
    void setLabel(const String& label)
    {
        label_ = label;
    }

    /// Returns the centroid m/z.
    DoubleReal getCentroidMZ()
    {
        return centroid_mz_;
    }

    DoubleReal getCentroidMZ() const
    {
        return centroid_mz_;
    }

    /// Returns the centroid RT.
    DoubleReal getCentroidRT()
    {
        return centroid_rt_;
    }

    DoubleReal getCentroidRT() const
    {
        return centroid_rt_;
    }

    DoubleReal getCentroidSD()
    {
        return centroid_sd_;
    }

    DoubleReal getCentroidSD() const
    {
        return centroid_sd_;
    }

    void setCentroidSD(const DoubleReal& tmp_sd)
    {
        centroid_sd_ = tmp_sd;
    }


    /// Gets smoothed intensities (empty if no smoothing was explicitly done beforehand!).
    std::vector<DoubleReal> getSmoothedIntensities()
    {
        return smoothed_intensities_;
    }

    std::vector<DoubleReal> getSmoothedIntensities() const
    {
        return smoothed_intensities_;
    }

    /// Set smoothed intensities (smoothing is done externally, e.g. by LowessSmoothing).
    void setSmoothedIntensities(const std::vector<DoubleReal>& db_vec)
    {
        if (trace_peaks_.size() != db_vec.size())
        {
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Number of smoothed intensities deviates from mass trace size! Aborting...", String(db_vec.size()));
        }

        smoothed_intensities_ = db_vec;
    }

    DoubleReal computeSmoothedPeakArea();
    DoubleReal computeSNR(bool, DoubleReal);

    /// Return estimated number of peaks spanning the full-width-at-half-maximum (previous estimation needed!).
    Size getFWHMScansNum()
    {
        return fwhm_num_scans_;
    }

    /// Set estimated number of peaks spanning the full-width-at-half-maximum.
    void setFWHMScansNum(Size r_fwhm)
    {
        fwhm_num_scans_ = r_fwhm;
    }

    /// Get scan time of mass trace
    DoubleReal getScanTime()
    {
        return scan_time_;
    }

    /** @name Computational methods
        */
    /// Sum up mass trace peak intensities for chromatographic peak area estimation.
    DoubleReal computePeakArea();
    DoubleReal computePeakArea(bool);

    /// Return the index of the mass trace's highest peak within the MassTrace container (based either on raw or smoothed intensities).
    Size findMaxByIntPeak(bool) const;

    /// Estimate FWHM of chromatographic peak in seconds (based on either raw or smoothed intensities). As a side-effect, the rough estimation of the number of scans within the FWHM range will be updated (see setFWHMScansNum).
    DoubleReal estimateFWHM(bool) const;

    /// Find local extrema within mass trace and return their indices.
    void findLocalExtrema(const Size&, std::vector<Size>&, std::vector<Size>&);

    /// Return the mass trace's convex hull.
    ConvexHull2D getConvexhull() const;


    /** @name Update methods for centroid RT and m/z
        */

    void updateSmoothedMaxRT();

    /// Compute & update centroid RT as a intensity-weighted mean of RTs.
    void updateWeightedMeanRT();

    /// Compute & update centroid RT as median position of intensities.
    void updateMedianRT();

    /// Compute & update centroid m/z as median of m/z values.
    void updateMedianMZ();

    /// Compute & update centroid m/z as mean of m/z values.
    void updateMeanMZ();

    /// Compute & update centroid m/z as weighted mean of m/z values.
    void updateWeightedMeanMZ();


private:
    /// Actual MassTrace container for doing centroid calculation, peak width estimation etc.
    std::vector<PeakType> trace_peaks_;

    /// Centroid m/z
    DoubleReal centroid_mz_;

    /// intensity-weighted STD
    DoubleReal centroid_sd_;

    /// Centroid RT
    DoubleReal centroid_rt_;

    /// Trace label
    String label_;

    /// Container for smoothed intensities. Smoothing must be done externally.
    std::vector<DoubleReal> smoothed_intensities_;

    /// Scan time (time difference between two consecutive scans)
    DoubleReal scan_time_;

    /// Rough estimate of a chromatographic peak's width (number of scans within the FWHM range).
    Size fwhm_num_scans_;
};

}

#endif // OPENMS_KERNEL_MASSTRACE_H

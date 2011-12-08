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

    /** @brief A container type that gathers peaks similar in m/z and moving along retention time. Depending on the method of extraction a mass trace could virtually represent
      a complete ion chromatogram (XIC) or merely a part of it (e.g., a chromatographic peak). The kernel class provides methods for computing mass trace characteristics such
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
        MassTrace(const std::list<PeakType>& );

        /// Detailed constructor 2
        MassTrace(const std::vector<PeakType>& );


        /// Default destructor
        ~MassTrace();
        /// Copy constructor
        MassTrace(const MassTrace &tr_obj);

        MassTrace& operator= (const MassTrace& rhs);

        /** @name Iterators
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

        /** @name Getters and Setters
        */

        /// Returns the number of peaks contained in the mass trace
        inline Size getSize() const
        {
            return trace_peaks_.size();
        }

        /// Get label of mass trace
        inline String getLabel()
        {
            return label_;
        }

        /// Set label of mass trace
        inline void setLabel(const String& label)
        {
            label_ = label;
        }

        /// Return the centroid m/z
        inline DoubleReal getCentroidMZ()
        {
            return centroid_mz_;
        }

        /// Return the centroid RT
        inline DoubleReal getCentroidRT()
        {
            return centroid_rt_;
        }

        /// Get smoothed intensities (empty if no smoothing was explicitly done beforehand!)
        inline std::vector<DoubleReal> getSmoothedIntensities()
        {
            return smoothed_intensities_;
        }

        inline std::vector<DoubleReal> getSmoothedIntensities() const
        {
            return smoothed_intensities_;
        }

        /// Set smoothed intensities
        inline void setSmoothedIntensities(const std::vector<DoubleReal>& db_vec)
        {
            if (trace_peaks_.size() != db_vec.size())
            {
                throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Number of smoothed intensities deviates from mass trace size! Aborting...", String(db_vec.size()));
            }

            smoothed_intensities_ = db_vec;
        }

        /// Return estimated number of peaks spanning the full-width-at-half-maximum (previous estimation needed!)
        inline Size getRoughFWHMsize()
        {
            return rough_fwhm_points_;
        }

        /// Set estimated number of peaks spanning the full-width-at-half-maximum
        inline void setRoughFWHMsize(Size r_fwhm)
        {
            rough_fwhm_points_ = r_fwhm;
        }




        /** @name Helper functions
        */
        /// Sum up mass trace peak intensities for chromatographic peak area estimation
        DoubleReal computePeakArea();

        /// Return a pointer to the mass trace's highest peak (based either on raw or smoothed intensities)
        Size findMaxByIntPeak(bool);

        /// Estimate FWHM of chromatographic peak (based on smoothed intensities)
        DoubleReal estimateFWHM();

        /// Find local extrema in mass trace
        void findLocalExtrema(const Size&, std::vector<Size>&, std::vector<Size>&);

        /// Return the mass trace's convex hull
        ConvexHull2D getConvexhull() const;


        /** @name Update methods for centroid RT and m/z
        */

        /// Compute & update centroid RT as a intensity-weighted mean of RTs
        void updateWeightedMeanRT();

        /// Compute & update centroid RT as median position of intensities
        void updateMedianRT();

        /// Compute & update centroid m/z as median of m/z values
        void updateMedianMZ();

        /// Compute & update centroid m/z as mean of m/z values
        void updateMeanMZ();

        /// Compute & update centroid m/z as weighted mean of m/z values
        void updateWeightedMeanMZ();


    private:
        /// Actual MassTrace container for doing centroid calculation, peak width estimation etc.
        std::vector<PeakType> trace_peaks_;

        /// Centroid m/z
        DoubleReal centroid_mz_;

        /// Centroid RT
        DoubleReal centroid_rt_;

        /// Trace label
        String label_;

        /// Container for smoothed intensities. Smoothing must be done externally.
        std::vector<DoubleReal> smoothed_intensities_;

        /// Rough estimate of a chromatographic peak's width (can be set while collecting the mass trace)
        Size rough_fwhm_points_;
    };

}

#endif // OPENMS_KERNEL_MASSTRACE_H

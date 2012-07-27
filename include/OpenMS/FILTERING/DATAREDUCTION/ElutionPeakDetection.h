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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#ifndef OPENMS_ElutionPeakDetection_H
#define OPENMS_ElutionPeakDetection_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>


namespace OpenMS
{
/**
  @brief Extracts chromatographic peaks from a mass trace.

    Mass traces may consist of several consequitively (partly overlapping) eluting peaks, e.g., stemming
    from isomeric compounds with exactly the same mass but different retentional behaviour. This method
    first applies LOWESS smoothing on the mass trace's intensities, then detects local minima/maxima in
    order to separate the chromatographic peaks from each other. This results in a vector that gathers
    the splitted mass traces (see @ref ElutionPeakDetection parameters).

  @htmlinclude OpenMS_ElutionPeakDetection.parameters

  @ingroup Quantitation
*/
    class OPENMS_DLLAPI ElutionPeakDetection
        : public DefaultParamHandler, public ProgressLogger
    {
    public:
        /// Default Constructor
        ElutionPeakDetection();

        /// Destructor
        virtual ~ElutionPeakDetection();

        /** @name Main computation methods
            */
        /// Extracts chromatographic peaks from a single MassTrace and stores the splits into a vector of new mass traces.
        void detectPeaks(MassTrace&, std::vector<MassTrace>&);

        /// Applies the aforementioned detection method on a series of mass traces as input.
        void detectPeaks(std::vector<MassTrace>&, std::vector<MassTrace>&);

        /// Computes an estimate of the average peak width of the experiment based on smoothed intensities (median) and an estimate of a lower and upper bound for the peak width (+/-2*MAD, median of absolute deviances).
        void filterByPeakWidth(std::vector<MassTrace>&, std::vector<MassTrace>&);


    protected:
        virtual void updateMembers_();

    private:
        DoubleReal chrom_fwhm_;
        // Size window_size_;
        DoubleReal sample_rate_;

        bool pw_filtering_;

        void detectElutionPeaks_(MassTrace&, std::vector<MassTrace>&);
    };
} // namespace OpenMS
#endif // OPENMS_ElutionPeakDetection_H

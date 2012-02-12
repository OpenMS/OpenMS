// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FILTERING_DATAREDUCTION_MASSTRACEDETECTION_H
#define OPENMS_FILTERING_DATAREDUCTION_MASSTRACEDETECTION_H

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
/**
  @brief A mass trace extraction method that gathers peaks similar in m/z and moving along retention time.

  Peaks of a @ref MSExperiment are sorted by their intensity and stored in a list of potential chromatographic
  apex positions. Starting with these, mass traces are extended in- and decreasingly in retention
  time. During this extension phase, the centroid m/z is computed on-line as an intensity-weighted mean of
  peaks. The extension phase ends when the frequency of gathered peaks drops below a
  threshold (min_sample_rate, see @ref MassTraceDetection parameters).

  @htmlinclude OpenMS_MassTraceDetection.parameters

  @ingroup Quantitation
*/


    class OPENMS_DLLAPI MassTraceDetection :
            public DefaultParamHandler,
            public ProgressLogger
    {
    public:
        /// Default constructor
        MassTraceDetection();

        /// Default destructor
        virtual ~MassTraceDetection();

    /** @name Helper methods
        */
    /// Allows the iterative computation of the intensity-weighted mean of a mass trace's centroid m/z.
    void updateIterativeWeightedMeanMZ(const DoubleReal&, const DoubleReal&, DoubleReal&, DoubleReal&, DoubleReal&);

    /// Computes a rough estimate of the average peak width of the experiment (median) and an estimate of a lower and upper bound for the peak width (+/-2*MAD, median of absolute deviances).
        void filterByPeakWidth(std::vector<MassTrace>&, std::vector<MassTrace>&);

    /** @name Main computation methods
        */
    /// Main method of MassTraceDetection. Extracts mass traces of a @ref MSExperiment and gathers them into a vector container.
        void run(const MSExperiment<Peak1D>&, std::vector<MassTrace>&);

    /// Invokes the run method (see above) on merely a subregion of a @ref MSExperiment map.
    void run(MSExperiment<Peak1D>::ConstAreaIterator& begin, MSExperiment<Peak1D>::ConstAreaIterator& end, std::vector<MassTrace>& found_masstraces);

    protected:
        virtual void updateMembers_();

    private:
    // parameter stuff
        DoubleReal mass_error_ppm_;
        DoubleReal noise_threshold_int_;
        DoubleReal chrom_apex_snt_;
        DoubleReal chrom_fwhm_;
        DoubleReal min_sample_rate_;
    };
}

#endif // OPENMS_FILTERING_DATAREDUCTION_MASSTRACEDETECTION_H

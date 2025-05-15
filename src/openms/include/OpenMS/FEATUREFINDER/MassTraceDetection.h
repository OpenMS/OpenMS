// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{

    /**
      @brief A mass trace extraction method that gathers peaks similar in m/z and moving along retention time.

      Peaks of a @ref MSExperiment are sorted by their intensity and stored in a
      list of potential chromatographic apex positions. Only peaks that are above
      the noise threshold (user-defined) are analyzed and only peaks that are n
      times above this minimal threshold are considered as apices. This saves
      computational resources and decreases the noise in the resulting output.

      Starting with these, mass traces are extended in- and decreasingly in
      retention time. During this extension phase, the centroid m/z is computed
      on-line as an intensity-weighted mean of peaks.

      The extension phase ends when either the frequency of gathered peaks drops
      below a threshold (min_sample_rate, see @ref MassTraceDetection parameters)
      or when the number of missed scans exceeds a threshold
      (trace_termination_outliers, see @ref MassTraceDetection parameters).

      Finally, only mass traces that pass a filter (a certain minimal and maximal
      length as well as having the minimal sample rate criterion fulfilled) get
      added to the result.

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
        ~MassTraceDetection() override;

        /** @name Helper methods
        */

        /// Allows the iterative computation of the intensity-weighted mean of a mass trace's centroid m/z.
        void updateIterativeWeightedMeanMZ(const double &, const double &, double &, double &, double &);

        /** @name Main computation methods
        */

        /// Main method of MassTraceDetection. Extracts mass traces of a @ref MSExperiment and gathers them into a vector container.
        void run(const PeakMap &, std::vector<MassTrace> &, const Size max_traces = 0);

        /// Invokes the run method (see above) on merely a subregion of a @ref MSExperiment map.
        void run(PeakMap::ConstAreaIterator & begin, PeakMap::ConstAreaIterator & end, std::vector<MassTrace> & found_masstraces);

        /** @name Private methods and members
        */
    protected:
        void updateMembers_() override;

    private:

        struct Apex
        {
          Apex(double intensity, Size scan_idx, Size peak_idx);
          double intensity;
          Size scan_idx;
          Size peak_idx;
        };

        /// The internal run method
        void run_(const std::vector<Apex>& chrom_apices,
                  const Size peak_count,
                  const PeakMap & work_exp,
                  const std::vector<Size>& spec_offsets,
                  std::vector<MassTrace> & found_masstraces,
                  const Size max_traces = 0);

        // parameter stuff
        double mass_error_ppm_;
        double noise_threshold_int_;
        double chrom_peak_snr_;
        MassTrace::MT_QUANTMETHOD quant_method_;

        String trace_termination_criterion_;
        Size trace_termination_outliers_;
        double min_sample_rate_;
        double min_trace_length_;
        double max_trace_length_;

        bool reestimate_mt_sd_;
    };
}

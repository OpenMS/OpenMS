// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken, Tristan Aretz, Manuel Zschaebitz $
// --------------------------------------------------------------------------

#pragma once

#include <boost/dynamic_bitset.hpp>
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
        void updateIterativeWeightedMeanMZ(const double added_mz,
                                           const double added_int,
                                           double& centroid_mz,
                                           double& prev_counter,
                                           double& prev_denom
                                          );

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
          Apex(PeakMap& map, const Size scan_idx, const Size peak_idx);
          // Default constructors
          Apex() = default;
          Apex(const Apex& other) = default;
          Apex(Apex&& other) = default;

          // Move assignment operator
          Apex& operator=(Apex&& other) = default;

          std::reference_wrapper<PeakMap> map_;
          Size scan_idx_;
          Size peak_idx_;

          ///get's the corresponding values
          double getMZ() const;
          double getRT() const;
          double getIntensity() const;
        };

        struct NextIndex
        {
          /// C'tor: init with number of threads in parallel region
          NextIndex(const std::vector<Apex>& data, const Size total_peak_count, const std::vector<Size>& spec_offsets, const double mass_error_ppm);

          /// Get the next free apex index which is not in the neighbourhood of a currently processing apex (in another thread)
          /// (Internally adds the apex's m/z to a blacklist which prevents other threads from obtaining an apex nearby)
          /// This function blocks until the next free apex is not conflicting anymore - i.e. another thread called setApexAsProcessed()
          Size getNextFreeIndex();
          
          /// If an apex was processed call this function to remove the apex from the blacklist and increase the current_apex_
          /// ... doesn't create a feature
          void setApexAsProcessed();
          /// ... does create a feature
          void setApexAsProcessed(const std::vector<std::pair<Size, Size> >& gathered_idx);

          bool isConflictingApex(const Apex a) const;

          bool isVisited(const Size scan_idx, const Size peak_idx) const;

          void setNumberOfThreads(const Size thread_num);


          /// reference for usage
          const std::vector<Apex>& data_;
          const std::vector<Size>& spec_offsets_;
        
          /// own datastructure
          std::vector<bool> peak_visited_;
          Size current_Apex_;
          std::vector<double> lock_list_;
          double mass_error_ppm_;
        };

        /// internal check for FWHM meta data
        bool checkFWHMMetaData_(const PeakMap& work_exp);

        /// The internal run method
        void run_(std::vector<Apex>& chrom_apices,
                  const Size peak_count,
                  const PeakMap & work_exp,
                  const std::vector<Size>& spec_offsets,
                  std::vector<MassTrace> & found_masstraces,
                  const Size max_traces = 0);

        // Find Offset for Peak
        static double findOffset_(const double centroid_mz, const double mass_error_ppm_);
        
        // parameter stuff
        double mass_error_ppm_;
        double mass_error_da_;
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
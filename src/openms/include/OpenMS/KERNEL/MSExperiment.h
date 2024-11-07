// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/AreaIterator.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>

#include <vector>


namespace OpenMS
{
  class Peak1D;
  class ChromatogramPeak;

  /**
    @brief In-Memory representation of a mass spectrometry run.

    This representation of an MS run is organized as list
    of spectra and chromatograms and provides an in-memory representation of
    popular mass-spectrometric file formats such as mzXML or mzML. The
    meta-data associated with an experiment is contained in
    ExperimentalSettings (by inheritance) while the raw data (as well as
    spectra and chromatogram level meta data) is stored in objects of type
    MSSpectrum and MSChromatogram, which are accessible through the getSpectrum
    and getChromatogram functions.

    @note For range operations, see \ref RangeUtils "RangeUtils module"!
    @note Some of the meta data is associated with the spectra directly (e.g. DataProcessing) and therefore the spectra need to be present to retain this information.
    @note For an on-disc representation of an MS experiment, see OnDiskExperiment.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI MSExperiment final : public RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity, RangeMobility>,
    public ExperimentalSettings
  {

public:
    typedef Peak1D PeakT;
    typedef ChromatogramPeak ChromatogramPeakT;

    /// @name Base type definitions
    //@{
    /// Peak type
    typedef PeakT PeakType;
    /// Chromatogram peak type
    typedef ChromatogramPeakT ChromatogramPeakType;
    /// Coordinate type of peak positions
    typedef PeakType::CoordinateType CoordinateType;
    /// Intensity type of peaks
    typedef PeakType::IntensityType IntensityType;
    /// RangeManager type
    typedef RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility> RangeManagerType;
    /// RangeManager type
    typedef RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity, RangeMobility> RangeManagerContainerType;
    /// Spectrum Type
    typedef MSSpectrum SpectrumType;
    /// Chromatogram type
    typedef MSChromatogram ChromatogramType;
    /// STL base class type
    typedef std::vector<SpectrumType> Base;
    //@}

    /// @name Iterator type definitions
    //@{
    /// Mutable iterator
    typedef std::vector<SpectrumType>::iterator Iterator;
    /// Non-mutable iterator
    typedef std::vector<SpectrumType>::const_iterator ConstIterator;
    /// Mutable area iterator type (for traversal of a rectangular subset of the peaks)
    typedef Internal::AreaIterator<PeakT, PeakT&, PeakT*, Iterator, SpectrumType::Iterator> AreaIterator;
    /// Immutable area iterator type (for traversal of a rectangular subset of the peaks)
    typedef Internal::AreaIterator<const PeakT, const PeakT&, const PeakT*, ConstIterator, SpectrumType::ConstIterator> ConstAreaIterator;
    //@}

    /// @name Delegations of calls to the vector of MSSpectra
    // Attention: these refer to the spectra vector only!
    //@{
    typedef Base::value_type value_type;
    typedef Base::iterator iterator;
    typedef Base::const_iterator const_iterator;

    /// Constructor
    MSExperiment();

    /// Copy constructor
    MSExperiment(const MSExperiment & source);

    /// Move constructor
    MSExperiment(MSExperiment&&) = default;

    /// Assignment operator
    MSExperiment & operator=(const MSExperiment & source);

    /// Move assignment operator
    MSExperiment& operator=(MSExperiment&&) & = default;

    /// Assignment operator
    MSExperiment & operator=(const ExperimentalSettings & source);

    /// D'tor
    ~MSExperiment() override;

    /// Equality operator
    bool operator==(const MSExperiment & rhs) const;

    /// Equality operator
    bool operator!=(const MSExperiment & rhs) const;
    
    /// The number of spectra
    inline Size size() const noexcept
    {
      return spectra_.size();
    }

    /// Resize to @p n spectra
    inline void resize(Size n)
    {
      spectra_.resize(n);
    }

    /// Are there any spectra (does not consider chromatograms)
    inline bool empty() const noexcept
    {
      return spectra_.empty();
    }
    
    /// Reserve space for @p n spectra
    inline void reserve(Size n)
    {
      spectra_.reserve(n);
    }

    /// Random access to @p n'th spectrum
    inline SpectrumType& operator[](Size n)
    {
      return spectra_[n];
    }

    /// Random access to @p n'th spectrum
    inline const SpectrumType& operator[](Size n) const
    {
      return spectra_[n];
    }

    inline Iterator begin() noexcept
    {
      return spectra_.begin();
    }

    inline ConstIterator begin() const noexcept
    {
      return spectra_.cbegin();
    }

    inline ConstIterator cbegin() const noexcept
    {
      return spectra_.cbegin();
    }

    inline Iterator end()
    {
      return spectra_.end();
    }

    inline ConstIterator end() const noexcept
    {
      return spectra_.cend();
    }
    
    inline ConstIterator cend() const noexcept
    {
      return spectra_.cend();
    }
    //@}

    // Aliases / chromatograms
    void reserveSpaceSpectra(Size s);
    void reserveSpaceChromatograms(Size s);

    ///@name Conversion to/from 2D data
    //@{
    /**
      @brief Reads out a 2D Spectrum

      Container can be a PeakArray or an STL container of peaks which
      supports push_back(), end() and back()
    */
    template <class Container>
    void get2DData(Container& cont) const
    {
      for (typename Base::const_iterator spec = spectra_.begin(); spec != spectra_.end(); ++spec)
      {
        if (spec->getMSLevel() != 1)
        {
          continue;
        }
        typename Container::value_type s; // explicit object here, since instantiation within push_back() fails on VS<12
        for (typename SpectrumType::const_iterator it = spec->begin(); it != spec->end(); ++it)
        {
          cont.push_back(s);
          cont.back().setRT(spec->getRT());
          cont.back().setMZ(it->getMZ());
          cont.back().setIntensity(it->getIntensity());
        }
      }
    }

    /**
      @brief Assignment of a data container with RT and MZ to an MSExperiment

      Fill MSExperiment with data.
      Note that all data present (including meta-data) will be deleted prior to adding new data!

      @param container An iterable type whose elements support getRT(), getMZ() and getIntensity()

      @exception Exception::Precondition is thrown if the container is not sorted according to
      retention time (in debug AND release mode)
    */
    template <class Container>
    void set2DData(const Container& container)
    {
      set2DData<false, Container>(container);
    }

    /**
      @brief Assignment of a data container with RT and MZ to an MSExperiment

      Fill MSExperiment with data.
      Note that all data present (including meta-data) will be deleted prior to adding new data!

      @param container An iterable type whose elements support getRT(), getMZ() and getIntensity()
      @param store_metadata_names [MetaInfoInterface input only] Names of metadata arrays which should be created;
                                  data is filled from the metainfointerface of each element of the input container.
                                  Currently, only float data is supported!

      @exception Exception::Precondition is thrown if the container is not sorted according to
      retention time (in debug AND release mode)
    */
    template <class Container>
    void set2DData(const Container& container, const StringList& store_metadata_names)
    {
      // clean up the container first
      clear(true);
      SpectrumType* spectrum = nullptr;
      typename PeakType::CoordinateType current_rt = -std::numeric_limits<typename PeakType::CoordinateType>::max();
      for (typename Container::const_iterator iter = container.begin(); iter != container.end(); ++iter)
      {
        // check if the retention time has changed
        if (current_rt != iter->getRT() || spectrum == nullptr)
        {
          // append new spectrum
          if (current_rt > iter->getRT())
          {
            throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input container is not sorted!");
          }
          current_rt =  iter->getRT();
          spectrum = createSpec_(current_rt, store_metadata_names);
        }

        // add either data point or mass traces (depending on template argument value)
        ContainerAdd_<typename Container::value_type, false>::addData_(spectrum, &(*iter), store_metadata_names);
      }
    }

     /**
      @brief Assignment of a data container with RT and MZ to an MSExperiment

      Fill MSExperiment with data.
      Note that all data present (including meta-data) will be deleted prior to adding new data!

      @tparam Container An iterable type whose elements support getRT(), getMZ() and getIntensity()
      @tparam add_mass_traces If true, each container element is searched for the metavalue
                             "num_of_masstraces".
                             If found, "masstrace_intensity" (X>=0) meta values are added as data points (with 13C spacing).
                             This is useful for, e.g., FF-Metabo output.
                             Note that the actual feature will NOT be added if mass traces are found (since MT0 is usually identical)
      @param container The input data with RT,m/z and intensity

      @exception Exception::Precondition is thrown if the container is not sorted according to
      retention time (in debug AND release mode) OR a "masstrace_intensity" value is expected but not found
         
    */
    template <bool add_mass_traces, class Container>
    void set2DData(const Container& container)
    {
      // clean up the container first
      clear(true);
      SpectrumType* spectrum = nullptr;
      typename PeakType::CoordinateType current_rt = -std::numeric_limits<typename PeakType::CoordinateType>::max();
      for (typename Container::const_iterator iter = container.begin(); iter != container.end(); ++iter)
      {
        // check if the retention time has changed
        if (current_rt != iter->getRT() || spectrum == nullptr)
        {
          // append new spectrum
          if (current_rt > iter->getRT())
          {
            throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input container is not sorted!");
          }
          current_rt =  iter->getRT();
          spectrum = createSpec_(current_rt);
        }

        // add either data point or mass traces (depending on template argument value)
        ContainerAdd_<typename Container::value_type, add_mass_traces>::addData_(spectrum, &(*iter));
      }
    }

    //@}


    ///@name Iterating ranges and areas
    //@{
    /// Returns an area iterator for @p area
    AreaIterator areaBegin(CoordinateType min_rt, CoordinateType max_rt, 
      CoordinateType min_mz, CoordinateType max_mz, UInt ms_level = 1);

    /// Returns an area iterator for all peaks in @p range. If a dimension is empty(), it is ignored (i.e. does not restrict the area)
    AreaIterator areaBegin(const RangeManagerType& range, UInt ms_level = 1);

    /// Returns an invalid area iterator marking the end of an area
    AreaIterator areaEnd();

    /// Returns a non-mutable area iterator for @p area
    ConstAreaIterator areaBeginConst(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, UInt ms_level = 1) const;

    /// Returns a non-mutable area iterator for all peaks in @p range. If a dimension is empty(), it is ignored (i.e. does not restrict the area)
    ConstAreaIterator areaBeginConst(const RangeManagerType& range, UInt ms_level = 1) const;

    /// Returns a non-mutable invalid area iterator marking the end of an area
    ConstAreaIterator areaEndConst() const;

    /* @brief Retrieves the peak data in the given mz-rt range and store data spectrum-wise in separate arrays.
     * 
     * For fast pyOpenMS access to peak data in format: [rt, [mz, intensity]]
     * 
     * @param min_rt The minimum retention time.
     * @param max_rt The maximum retention time.
     * @param min_mz The minimum m/z value.
     * @param max_mz The maximum m/z value.
     * @param ms_level The MS level of the spectra to consider.
     * @param rt The vector to store the retention times in.
     * @param mz The vector to store the m/z values in.
     * @param intensity The vector to store the intensities in.
     */
    void get2DPeakDataPerSpectrum(
      CoordinateType min_rt, 
      CoordinateType max_rt, 
      CoordinateType min_mz, 
      CoordinateType max_mz,
      Size ms_level,
      std::vector<float>& rt, 
      std::vector<std::vector<float>>& mz, 
      std::vector<std::vector<float>>& intensity) const
    {
      float t = -1.0;
      for (auto it = areaBeginConst(min_rt, max_rt, min_mz, max_mz, ms_level); it != areaEndConst(); ++it)
      {
        if (it.getRT() != t) 
        {
          t = (float)it.getRT();
          rt.push_back(t);
          mz.push_back(std::vector<float>());
          intensity.push_back(std::vector<float>());
        }
        mz.back().push_back((float)it->getMZ());
        intensity.back().push_back(it->getIntensity());
      }
    }

    /* @brief Retrieves the peak data in the given mz-rt range and store data spectrum-wise in separate arrays.
     * 
     * For fast pyOpenMS access to MS1 peak data in format: [rt, [mz, intensity, ion mobility]]
     * 
     * @param min_rt The minimum retention time.
     * @param max_rt The maximum retention time.
     * @param min_mz The minimum m/z value.
     * @param max_mz The maximum m/z value.
     * @param ms_level The MS level of the spectra to consider.
     * @param rt The vector to store the retention times in.
     * @param mz The vector to store the m/z values in.
     * @param intensity The vector to store the intensities in.
     * @param ion_mobility The vector to store the ion mobility values in.
    */
    void get2DPeakDataIMPerSpectrum(
      CoordinateType min_rt, 
      CoordinateType max_rt, 
      CoordinateType min_mz, 
      CoordinateType max_mz,
      Size ms_level,     
      std::vector<float>& rt, 
      std::vector<std::vector<float>>& mz,
      std::vector<std::vector<float>>& intensity, 
      std::vector<std::vector<float>>& ion_mobility) const
    {
      DriftTimeUnit unit;
      std::vector<float> im;
      float t = -1.0;
      for (auto it = areaBeginConst(min_rt, max_rt, min_mz, max_mz, ms_level); it != areaEndConst(); ++it)
      {
        if (it.getRT() != t)
        {
          t = (float)it.getRT();
          rt.push_back(t);
          std::tie(unit, im) = it.getSpectrum().maybeGetIMData();
          mz.push_back(std::vector<float>());
          intensity.push_back(std::vector<float>());
          ion_mobility.push_back(std::vector<float>());
        }

        if (unit != DriftTimeUnit::NONE)
        {
          const Size peak_index = it.getPeakIndex().peak;
          ion_mobility.back().push_back(im[peak_index]);
        }
        else
        {
          ion_mobility.back().push_back(-1.0);
        }
        mz.back().push_back((float)it->getMZ());
        intensity.back().push_back(it->getIntensity());
      }
    }

    /* @brief Retrieves the peak data in the given mz-rt range and store in separate arrays.
     * 
     * For fast pyOpenMS access to MS1 peak data in format: [rt, mz, intensity]
     * 
     * @param min_rt The minimum retention time.
     * @param max_rt The maximum retention time.
     * @param min_mz The minimum m/z value.
     * @param max_mz The maximum m/z value.
     * @param ms_level The MS level of the spectra to consider.
     * @param rt The vector to store the retention times in.
     * @param mz The vector to store the m/z values in.
     * @param intensity The vector to store the intensities in.
    */    
    void get2DPeakData(
      CoordinateType min_rt,
      CoordinateType max_rt,
      CoordinateType min_mz,
      CoordinateType max_mz,
      Size ms_level,
      std::vector<float>& rt,
      std::vector<float>& mz,
      std::vector<float>& intensity) 
      const
    {
      for (auto it = areaBeginConst(min_rt, max_rt, min_mz, max_mz, ms_level); it != areaEndConst(); ++it)
      {
        rt.push_back((float)it.getRT());
        mz.push_back((float)it->getMZ());
        intensity.push_back(it->getIntensity());
      }
    }


    /* @brief Retrieves the peak data in the given mz-rt range and store in separate arrays.
     * 
     * For fast pyOpenMS access to MS1 peak data in format: [rt, mz, intensity, ion mobility]
     * 
     * @param min_rt The minimum retention time.
     * @param max_rt The maximum retention time.
     * @param min_mz The minimum m/z value.
     * @param max_mz The maximum m/z value.
     * @param ms_level The MS level of the spectra to consider.
     * @param rt The vector to store the retention times in.
     * @param mz The vector to store the m/z values in.
     * @param intensity The vector to store the intensities in.
    */
    void get2DPeakDataIM(
      CoordinateType min_rt,
      CoordinateType max_rt,
      CoordinateType min_mz,
      CoordinateType max_mz,
      Size ms_level,
      std::vector<float>& rt,
      std::vector<float>& mz,
      std::vector<float>& intensity,
      std::vector<float>& ion_mobility) const
    {
      for (auto it = areaBeginConst(min_rt, max_rt, min_mz, max_mz, ms_level); it != areaEndConst(); ++it)
      {
        DriftTimeUnit unit;
        std::vector<float> im;
        float t = -1.0;
        if (it.getRT() != t)
        {
          t = (float)it.getRT();
          std::tie(unit, im) = it.getSpectrum().maybeGetIMData();
        }
        rt.push_back((float)it.getRT());
        mz.push_back((float)it->getMZ());
        intensity.push_back(it->getIntensity());
        if (unit != DriftTimeUnit::NONE)
        {
          const Size peak_index = it.getPeakIndex().peak;
          ion_mobility.push_back(im[peak_index]);
        }
        else
        {
          ion_mobility.push_back(-1.0);
        }
      }
    }

  /**
   * @brief Calculates the sum of intensities for a range of elements.
   * 
   * @tparam Iterator The iterator type.
   * @param begin The iterator pointing to the beginning of the range.
   * @param end The iterator pointing to the end of the range.
   * @return The sum of intensities.
   * 
   * @throws static assert fails if the iterator value type does not have a `getIntensity()` member function.
   */
struct SumIntensityReduction {

  template <typename Iterator>
  auto operator()(Iterator begin, Iterator end) const {
    // Static assert to verify iterator type has intensity accessor
    using ValueType = typename std::iterator_traits<Iterator>::value_type;
    using IntensityType = decltype(std::declval<ValueType>().getIntensity());
    static_assert(std::is_member_function_pointer_v<decltype(&ValueType::getIntensity)>,
           "Iterator value type must have getIntensity() member function");

    IntensityType sum{};
    for (auto it = begin; it != end; ++it) {
      sum += it->getIntensity();
    }
    return sum;
  }
};

/**
 * @brief Aggregates data over specified m/z and RT ranges at a given MS level using a custom reduction function.
 *
 * This function processes spectra at a specified MS level and aggregates data within specified m/z (mass-to-charge)
 * and RT (retention time) ranges. For each m/z and RT range, it computes a result using the provided m/z reduction
 * function over the peaks that fall within the range.
 *
 * The results are organized in a two-dimensional vector, where each sub-vector corresponds to a particular m/z and RT
 * range, and contains the aggregated results for each spectrum that falls within that RT range.
 * 
 * This function allows for flexible aggregation of data within specified m/z and RT ranges, and is useful for XIC extraction.
 *
 * @tparam MzReductionFunctionType
 *   A callable type (function, lambda, or functor) that takes two iterators (`begin_it` and `end_it`) over peaks
 *   (`MSSpectrum::ConstIterator`) and returns a `CoordinateType`. The function defines how to reduce or aggregate
 *   the peaks within the specified m/z range (e.g., summing intensities, computing the mean m/z, etc.).
 *
 * @param[in,out] mz_rt_ranges
 *   A vector of pairs of `RangeMZ` and `RangeRT` specifying the m/z and RT ranges over which to aggregate data.
 *   Each pair defines a rectangular region in the m/z-RT plane. The vector will be sorted in-place by ascending
 *   minimum m/z and descending maximum m/z within the function.
 *
 * @param[in] ms_level
 *   The MS level of the spectra to be processed. Only spectra matching this MS level will be considered in the aggregation.
 *
 * @param[in] func_mz_reduction
 *   A function or functor that performs the m/z reduction. It should have the signature:
 *   `CoordinateType func_mz_reduction(MSSpectrum::ConstIterator begin_it, MSSpectrum::ConstIterator end_it);`
 *   The function receives iterators to the beginning and end of the peaks within the m/z range for a specific spectrum.
 *   It should process these peaks and return a single `CoordinateType` value representing the aggregated result.
 *
 * @return
 *   A vector of vectors of `MSExperiment::CoordinateType`. Each sub-vector corresponds to an m/z and RT range in
 *   `mz_rt_ranges`, and contains the aggregated results for each spectrum that falls within that RT range.
 *   The size of each sub-vector equals the number of spectra that overlap with the RT range.
 *
 * @note
 * - If `mz_rt_ranges` is empty or there are no spectra at the specified MS level, the function returns an empty vector.
 * - The `mz_rt_ranges` vector will be sorted within the function by ascending minimum m/z and descending maximum m/z.
 * - The function uses OpenMP for parallelization over spectra. Ensure that your reduction function is thread-safe.
 * - The aggregation is performed only on the peaks that fall within both the specified m/z and RT ranges.
 * - This methods works best with larger number of m/z and RT ranges and a large number of spectra.
 *
 * @warning
 * - The function modifies `mz_rt_ranges` by sorting it. If the original order is important, make a copy before calling.
 * - The provided `func_mz_reduction` must be able to handle empty ranges (i.e., when `begin_it == end_it`).
 *
 * @exception None
 */
template<class MzReductionFunctionType>
std::vector<std::vector<MSExperiment::CoordinateType>> aggregate(
    const std::vector<std::pair<RangeMZ, RangeRT>>& mz_rt_ranges,
    unsigned int ms_level,
    MzReductionFunctionType func_mz_reduction) const
{
    // Early exit if there are no ranges
    if (mz_rt_ranges.empty()) 
    {
      // likely an error, but we return an empty vector instead of throwing an exception for now
      return {};
    }

    // Create a view of the spectra with given MS level
    std::vector<std::reference_wrapper<const MSSpectrum>> spectra_view;
    spectra_view.reserve(spectra_.size());
    std::copy_if(spectra_.begin(), spectra_.end(), 
                 std::back_inserter(spectra_view),
                 [ms_level](const auto& spec) { 
                     return spec.getMSLevel() == ms_level; 
                 });

    // Early exit if there are no spectra with the given MS level
    if (spectra_view.empty()) { // could be valid use or an error -> we return an empty vector
      return {};
    }

    // Get the indices of the spectra covered by the RT ranges by considering the MS level
    // If start and stop are the same, the range is empty
    auto getCoveredSpectra = [](
        const std::vector<std::reference_wrapper<const MSSpectrum>>& spectra_view, 
        const std::vector<std::pair<RangeMZ, RangeRT>>& mz_rt_ranges) 
        ->  std::vector<std::pair<size_t, size_t>>
      {
        std::vector<std::pair<size_t, size_t>> res;
        res.reserve(mz_rt_ranges.size());
        
        for (const auto & mz_rt : mz_rt_ranges) 
        {
          // std::cout << "rt range: " << mz_rt.second.getMin() << " - " << mz_rt.second.getMax() << std::endl;
          // std::cout << "specs start:" << spectra_view[0].get().getRT() << " specs end:" << spectra_view[spectra_view.size() - 1].get().getRT() << std::endl;
          auto start_it = std::lower_bound(spectra_view.begin(), spectra_view.end(), mz_rt.second.getMin(), 
            [](const auto& spec, double rt) 
            { return spec.get().getRT() < rt; });

          auto stop_it = std::upper_bound(spectra_view.begin(), spectra_view.end(), mz_rt.second.getMax(), 
            [](double rt, const auto& spec) 
            { return rt < spec.get().getRT(); });

            res.emplace_back(
                std::distance(spectra_view.begin(), start_it),
                std::distance(spectra_view.begin(), stop_it)
            );
            // std::cout << "start: " << std::distance(spectra_view.begin(), start_it) << " stop: " << std::distance(spectra_view.begin(), stop_it) << std::endl;
        }
        return res;
      };

    // For each range, gets (spectrum start index, spectrum stop index). The spectra covered by each RT range.
    const std::vector<std::pair<size_t, size_t>> rt_ranges_idcs = getCoveredSpectra(spectra_view, mz_rt_ranges);

    // Initialize result vector
    std::vector<std::vector<MSExperiment::CoordinateType>> result(mz_rt_ranges.size());

    // Initialize counts per spectrum index and total mappings
    std::vector<std::vector<size_t>> spec_idx_to_range_idx(spectra_view.size());

    // Build spectrum to range index mapping
    for (size_t i = 0; i < rt_ranges_idcs.size(); ++i) 
    {
        const auto& [start, stop] = rt_ranges_idcs[i];
        result[i].resize(stop - start);
        // std::cout << "start: " << start << " stop: " << stop << std::endl;
        for (size_t j = start; j < stop; ++j) 
        {
            spec_idx_to_range_idx[j].push_back(i);
        }
    }

   #pragma omp parallel for schedule(dynamic)
   for (Int64 i = 0; i < (Int64)spec_idx_to_range_idx.size(); ++i) // OpenMP on windows still requires signed loop variable
   {
      if (spec_idx_to_range_idx[i].empty()) continue; // no ranges for this spectrum? skip it

      const auto& spec = spectra_view[i].get();
      auto spec_begin = spec.cbegin();
      auto spec_end = spec.cend();

      for (size_t range_idx : spec_idx_to_range_idx[i]) 
      {
        const auto& mz_range = mz_rt_ranges[range_idx].first;
        
        // Find data points within MZ range
        auto start_it = spec.PosBegin(spec_begin, mz_range.getMinMZ(), spec_end);
        auto end_it = start_it;
        
        while (end_it != spec_end && end_it->getPosition() <= mz_range.getMaxMZ()) 
        {
          ++end_it;
        }

        // std::cout << "calculating reduction on range: " << range_idx << " for spectrum: " << i << " and peaks " << std::distance(spec.begin(), start_it) << " - " << std::distance(spec.begin(), end_it) << std::endl;

        // Calculate result using provided reduction function
        result[range_idx][i - rt_ranges_idcs[range_idx].first] = 
            func_mz_reduction(start_it, end_it);
      }
    }
    return result;
  }

// Overload without func_mz_reduction parameter (default to SumIntensityReduction). Needed because of template deduction issues
std::vector<std::vector<MSExperiment::CoordinateType>> aggregate(
    const std::vector<std::pair<RangeMZ, RangeRT>>& mz_rt_ranges,
    unsigned int ms_level) const
{
  return aggregate(mz_rt_ranges, ms_level, SumIntensityReduction());
}

/**
 * @brief Extracts extracted ion chromatograms (XICs) from the MSExperiment.
 *
 * This function takes a vector of mz_rt_ranges, an ms_level, and a MzReductionFunctionType
 * and extracts the XICs from the MSExperiment based on the given parameters.
 *
 * @param mz_rt_ranges A vector of pairs of RangeMZ and RangeRT representing the m/z and retention time ranges.
 * @param ms_level The MS level of the spectra to consider.
 * @param func_mz_reduction The MzReductionFunctionType used to reduce the m/z values.
 *
 * @return A vector of MSChromatogram objects representing the extracted XICs.
 */
template<class MzReductionFunctionType>
std::vector<MSChromatogram> extractXICs(
    const std::vector<std::pair<RangeMZ, RangeRT>>& mz_rt_ranges,
    unsigned int ms_level,
    MzReductionFunctionType func_mz_reduction) const
{
    // Early exit if there are no ranges
    if (mz_rt_ranges.empty()) 
    {
      // likely an error, but we return an empty vector instead of throwing an exception for now
      return {};
    }

    // Create a view of the spectra with given MS level
    std::vector<std::reference_wrapper<const MSSpectrum>> spectra_view;
    spectra_view.reserve(spectra_.size());
    std::copy_if(spectra_.begin(), spectra_.end(), 
                 std::back_inserter(spectra_view),
                 [ms_level](const auto& spec) { 
                     return spec.getMSLevel() == ms_level; 
                 });

    // Early exit if there are no spectra with the given MS level
    if (spectra_view.empty()) { // could be valid use or an error -> we return an empty vector
      return {};
    }

    // Get the indices of the spectra covered by the RT ranges by considering the MS level
    // If start and stop are the same, the range is empty
    auto getCoveredSpectra = [](
        const std::vector<std::reference_wrapper<const MSSpectrum>>& spectra_view, 
        const std::vector<std::pair<RangeMZ, RangeRT>>& mz_rt_ranges) 
        ->  std::vector<std::pair<size_t, size_t>>
      {
        std::vector<std::pair<size_t, size_t>> res;
        res.reserve(mz_rt_ranges.size());

        for (const auto & mz_rt : mz_rt_ranges) 
        {
          auto start_it = std::lower_bound(spectra_view.begin(), spectra_view.end(), mz_rt.second.getMin(), 
            [](const auto& spec, double rt) 
            { return spec.get().getRT() < rt; });

          auto stop_it = std::upper_bound(spectra_view.begin(), spectra_view.end(), mz_rt.second.getMax(), 
            [](double rt, const auto& spec) 
            { return rt < spec.get().getRT(); });

            res.emplace_back(
                std::distance(spectra_view.begin(), start_it),
                std::distance(spectra_view.begin(), stop_it)
            );
        }
        return res;
      };

    // For each range, gets (spectrum start index, spectrum stop index). The spectra covered by each RT range.
    const std::vector<std::pair<size_t, size_t>> rt_ranges_idcs = getCoveredSpectra(spectra_view, mz_rt_ranges);

    // Initialize result vector
    std::vector<MSChromatogram> result(mz_rt_ranges.size());

    // Initialize counts per spectrum index and total mappings
    std::vector<std::vector<size_t>> spec_idx_to_range_idx(spectra_view.size());

    // Build spectrum to range index mapping
    for (size_t i = 0; i < rt_ranges_idcs.size(); ++i) 
    {        
        const auto& [start, stop] = rt_ranges_idcs[i];
        result[i].resize(stop - start);
        result[i].getProduct().setMZ(
          (mz_rt_ranges[i].first.getMinMZ() + mz_rt_ranges[i].first.getMaxMZ()) / 2.0);
        for (size_t j = start; j < stop; ++j) 
        {
          spec_idx_to_range_idx[j].push_back(i);
        }
    }

   #pragma omp parallel for schedule(dynamic)
   for (Int64 i = 0; i < (Int64)spec_idx_to_range_idx.size(); ++i) // OpenMP on windows still requires signed loop variable
   {
      if (spec_idx_to_range_idx[i].empty()) continue; // no ranges for this spectrum? skip 
      
      const auto& spec = spectra_view[i].get();
      const double rt = spec.getRT();
      auto spec_begin = spec.cbegin();
      auto spec_end = spec.cend();

      for (size_t range_idx : spec_idx_to_range_idx[i]) 
      {
        const auto& mz_range = mz_rt_ranges[range_idx].first;
        
        // Find data points within MZ range
        auto start_it = spec.PosBegin(spec_begin, mz_range.getMinMZ(), spec_end);
        auto end_it = start_it;
        
        while (end_it != spec_end && end_it->getPosition() <= mz_range.getMaxMZ()) 
        {
          ++end_it;
        }

        // Calculate result using provided reduction function
        result[range_idx][i - rt_ranges_idcs[range_idx].first] = 
            ChromatogramPeak(rt, func_mz_reduction(start_it, end_it));
      }
    }

    for (auto& r : result) r.updateRanges(); // TODO: prob.. faster to look at first and last peaks as range is sorted

    return result;
  }

// Overload without func_mz_reduction parameter (needed because of template deduction issue)
std::vector<MSChromatogram> extractXICs(
    const std::vector<std::pair<RangeMZ, RangeRT>>& mz_rt_ranges,
    unsigned int ms_level) const
{
    return extractXICs(mz_rt_ranges, ms_level, SumIntensityReduction());
}

  /**
   * @brief Wrapper for aggregate function that takes a matrix of m/z and RT ranges
   * 
   * @param ranges Matrix where each row contains [mz_min, mz_max, rt_min, rt_max]
   * @param ms_level MS level to process
   * @param mz_agg Aggregation function for m/z values ("sum", "max", "min", "mean")
   * @return Vector of vectors containing aggregated intensity values for each range
   */
  std::vector<std::vector<MSExperiment::CoordinateType>> aggregateFromMatrix(
      const Matrix<double>& ranges,
      unsigned int ms_level,
      const std::string& mz_agg) const
  {
      // Check matrix dimensions
      if (ranges.cols() != 4)
      {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Range matrix must have 4 columns [mz_min, mz_max, rt_min, rt_max]");
      }

      // Convert matrix rows to vector of pairs
      std::vector<std::pair<RangeMZ, RangeRT>> mz_rt_ranges;
      mz_rt_ranges.reserve((Size)ranges.rows());
      
      for (Size i = 0; i < (Size)ranges.rows(); ++i)
      {
          mz_rt_ranges.emplace_back(
              RangeMZ(ranges(i, 0), ranges(i, 1)), // min max mz
              RangeRT(ranges(i, 2), ranges(i, 3))  // min max rt
          );
          // std::cout << "mz: " << ranges(i, 0) << " - " << ranges(i, 1) << " rt: " << ranges(i, 2) << " - " << ranges(i, 3) << std::endl;
      }

      // Call appropriate aggregation function based on mz_agg parameter
      if (mz_agg == "sum")
      {
          return aggregate(mz_rt_ranges, ms_level,
              [](auto begin_it, auto end_it)
              {
                  return std::accumulate(begin_it, end_it, 0.0,
                      [](double a, const Peak1D& b) { return a + b.getIntensity(); });
              });
      }
      else if (mz_agg == "max")
      {
          return aggregate(mz_rt_ranges, ms_level,
              [](auto begin_it, auto end_it)->double
              {
                  if (begin_it == end_it) return 0.0;
                  return std::max_element(begin_it, end_it,
                      [](const Peak1D& a, const Peak1D& b) { return a.getIntensity() < b.getIntensity(); }
                  )->getIntensity();
              });
      }
      else if (mz_agg == "min")
      {
          return aggregate(mz_rt_ranges, ms_level,
              [](auto begin_it, auto end_it)->double
              {
                  if (begin_it == end_it) return 0.0;
                  return std::min_element(begin_it, end_it,
                      [](const Peak1D& a, const Peak1D& b) { return a.getIntensity() < b.getIntensity(); }
                  )->getIntensity();
              });
      }
      else if (mz_agg == "mean")
      {
          return aggregate(mz_rt_ranges, ms_level,
              [](auto begin_it, auto end_it)
              {
                  if (begin_it == end_it) return 0.0;
                  double sum = std::accumulate(begin_it, end_it, 0.0,
                      [](double a, const Peak1D& b) { return a + b.getIntensity(); });
                  return sum / static_cast<double>(std::distance(begin_it, end_it));
              });
      }
      else
      {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Invalid aggregation function", mz_agg);
      }
  }

  /**
   * @brief Wrapper for extractXICs function that takes a matrix of m/z and RT ranges
   * 
   * @param ranges Matrix where each row contains [mz_min, mz_max, rt_min, rt_max]
   * @param ms_level MS level to process
   * @param mz_agg Aggregation function for m/z values ("sum", "max", "min", "mean")
   * @return Vector of MSChromatogram objects, one for each range
   */
  std::vector<MSChromatogram> extractXICsFromMatrix(
      const Matrix<double>& ranges,
      unsigned int ms_level,
      const std::string& mz_agg) const
  {
      // Check matrix dimensions
      if (ranges.cols() != 4)
      {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Range matrix must have 4 columns [mz_min, mz_max, rt_min, rt_max]");
      }

      // Convert matrix rows to vector of pairs
      std::vector<std::pair<RangeMZ, RangeRT>> mz_rt_ranges;
      mz_rt_ranges.reserve((Size)ranges.rows());
      
      for (Size i = 0; i < (Size)ranges.rows(); ++i)
      {
          mz_rt_ranges.emplace_back(
              RangeMZ(ranges(i, 0), ranges(i, 1)),
              RangeRT(ranges(i, 2), ranges(i, 3))
          );
      }

      // Call appropriate extractXICs function based on mz_agg parameter
      if (mz_agg == "sum")
      {
          return extractXICs(mz_rt_ranges, ms_level,
              [](auto begin_it, auto end_it)
              {
                  return std::accumulate(begin_it, end_it, 0.0,
                      [](double a, const Peak1D& b) { return a + b.getIntensity(); });
              });
      }
      else if (mz_agg == "max")
      {
          return extractXICs(mz_rt_ranges, ms_level,
              [](auto begin_it, auto end_it)->double
              {
                  if (begin_it == end_it) return 0.0;
                  return std::max_element(begin_it, end_it,
                      [](const Peak1D& a, const Peak1D& b) { return a.getIntensity() < b.getIntensity(); }
                  )->getIntensity();
              });
      }
      else if (mz_agg == "min")
      {
          return extractXICs(mz_rt_ranges, ms_level,
              [](auto begin_it, auto end_it)->double
              {
                  if (begin_it == end_it) return 0.0;
                  return std::min_element(begin_it, end_it,
                      [](const Peak1D& a, const Peak1D& b) { return a.getIntensity() < b.getIntensity(); }
                  )->getIntensity();
              });
      }
      else if (mz_agg == "mean")
      {
          return extractXICs(mz_rt_ranges, ms_level,
              [](auto begin_it, auto end_it)
              {
                  if (begin_it == end_it) return 0.0;
                  double sum = std::accumulate(begin_it, end_it, 0.0,
                      [](double a, const Peak1D& b) { return a + b.getIntensity(); });
                  return sum / static_cast<double>(std::distance(begin_it, end_it));
              });
      }
      else
      {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Invalid aggregation function", mz_agg);
      }
  }

    /**
      @brief Fast search for spectrum range begin

      Returns the first scan which has equal or higher (>=) RT than @p rt.

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    ConstIterator RTBegin(CoordinateType rt) const;

    /**
      @brief Fast search for spectrum range end (returns the past-the-end iterator)

      Returns the first scan which has higher (>) RT than @p rt.

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    ConstIterator RTEnd(CoordinateType rt) const;

    /**
      @brief Fast search for spectrum range begin

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    Iterator RTBegin(CoordinateType rt);

    /**
      @brief Fast search for spectrum range end (returns the past-the-end iterator)

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    Iterator RTEnd(CoordinateType rt);


    /**
      @brief Fast search for spectrum range begin

      Returns the first scan which has equal or higher (>=) ion mobility than @p rt.

      @note Make sure the spectra are sorted with respect to ion mobility! Otherwise the result is undefined.
    */
    ConstIterator IMBegin(CoordinateType im) const;

    /**
      @brief Fast search for spectrum range end (returns the past-the-end iterator)

      Returns the first scan which has higher (>) ion mobility than @p im.

      @note Make sure the spectra are sorted with respect to ion mobility! Otherwise the result is undefined.
    */
    ConstIterator IMEnd(CoordinateType im) const;
    //@}

    /**
      @name Range methods

      @note The range values (min, max, etc.) are not updated automatically. Call updateRanges() to update the values!
    */
    ///@{
    // Docu in base class
    void updateRanges() override;

    /**
      @brief Updates the m/z, intensity, retention time and MS level ranges of all spectra with a certain ms level

      @param ms_level MS level to consider for m/z range , RT range and intensity range (All MS levels if negative)
    */
    void updateRanges(Int ms_level);

    /// returns the total number of peaks
    UInt64 getSize() const;

    /// returns an array of MS levels
    const std::vector<UInt>& getMSLevels() const;

    ///@}

    /// If the file is loaded from an sqMass file, this run-ID allows to connect to the corresponding OSW identification file
    /// If the run-ID was not stored (older version) or this MSExperiment was not loaded from sqMass, then 0 is returned.
    UInt64 getSqlRunID() const;

    /// sets the run-ID which is used when storing an sqMass file
    void setSqlRunID(UInt64 id);

    ///@name Sorting spectra and peaks
    ///@{
    /**
      @brief Sorts the data points by retention time

      @param sort_mz if @em true, spectra are sorted by m/z position as well
    */
    void sortSpectra(bool sort_mz = true);

    /**
      @brief Sorts the data points of the chromatograms by m/z

      @param sort_rt if @em true, chromatograms are sorted by rt position as well
    */
    void sortChromatograms(bool sort_rt = true);

    /**
      @brief Checks if all spectra are sorted with respect to ascending RT

      @param check_mz if @em true, checks if all peaks are sorted with respect to ascending m/z
    */
    bool isSorted(bool check_mz = true) const;

    //@}

    /// Clear all internal data (spectra, ranges, metadata)
    void reset();

    /**
      @brief Clears the meta data arrays of all contained spectra (float, integer and string arrays)

      @return @em true if meta data arrays were present and removed. @em false otherwise.
    */
    bool clearMetaDataArrays();

    /// returns the meta information of this experiment (const access)
    const ExperimentalSettings& getExperimentalSettings() const;

    /// returns the meta information of this experiment (mutable access)
    ExperimentalSettings& getExperimentalSettings();

    /// get the file path to the first MS run
    void getPrimaryMSRunPath(StringList& toFill) const;

    /**
      @brief Returns the precursor spectrum of the scan pointed to by @p iterator

      If there is no (matching) precursor scan the past-the-end iterator is returned.

      This assumes that precursors occur somewhere before the current spectrum
      but not necessarily the first one from the last MS level (we double-check with
      the annotated precursorList.

      If precursor annotations are present, uses the native spectrum ID from the 
      @em first precursor entry of the current scan
      for comparisons -> Works for multiple precursor ranges from the same precursor scan
      but not for multiple precursor ranges from different precursor scans.
      If none are present, picks the first scan of a lower level.
    */
    ConstIterator getPrecursorSpectrum(ConstIterator iterator) const;

    /**
      @brief Returns the index of the precursor spectrum for spectrum at index @p zero_based_index

      If there is no precursor scan -1 is returned. 
    */
    int getPrecursorSpectrum(int zero_based_index) const;

    /**
      @brief Returns the first product spectrum of the scan pointed to by @p iterator
      A product spectrum is a spectrum of the next higher MS level that has the
      current spectrum as precursor.
      If there is no product scan, the past-the-end iterator is returned.
      This assumes that product occurs somewhere after the current spectrum
      and comes before the next scan that is of a level that is lower than
      the current one.
\verbatim
      Example:
      MS1 - ix: 0
        MS2 - ix: 1, prec: 0
        MS2 - ix: 2, prec: 0 <-- current scan
        MS3 - ix: 3, prec: 1
        MS3 - ix: 4, prec: 2 <-- product scan
        MS2 - ix: 5, prec: 0
        MS3 - ix: 6, prec: 5
      MS1 - ix: 7
        ...  <-- Not searched anymore. Returns end of experiment iterator if not found until here.
\endverbatim
      Uses the native spectrum ID from the @em first precursor entry of the potential product scans
      for comparisons -> Works for multiple precursor ranges from the same precursor scan
      but not for multiple precursor ranges from different precursor scans.
    */
    ConstIterator getFirstProductSpectrum(ConstIterator iterator) const;

    /**
      @brief Returns the index of the first product spectrum given an index.

      @param zero_based_index The index of the current spectrum.

      @return Index of the first product spectrum or -1 if not found.
    */
    int getFirstProductSpectrum(int zero_based_index) const;

    /// Swaps the content of this map with the content of @p from
    void swap(MSExperiment& from);

    /// sets the spectrum list
    void setSpectra(const std::vector<MSSpectrum>& spectra);
    void setSpectra(std::vector<MSSpectrum>&& spectra);

    /// adds a spectrum to the list
    void addSpectrum(const MSSpectrum& spectrum);
    void addSpectrum(MSSpectrum&& spectrum);

    /// returns the spectrum list
    const std::vector<MSSpectrum>& getSpectra() const;

    /// returns the spectrum list (mutable)
    std::vector<MSSpectrum>& getSpectra();

    /// Returns the closest(=nearest) spectrum in retention time to the given RT
    ConstIterator getClosestSpectrumInRT(const double RT) const;
    Iterator getClosestSpectrumInRT(const double RT);

    /// Returns the closest(=nearest) spectrum in retention time to the given RT of a certain MS level
    ConstIterator getClosestSpectrumInRT(const double RT, UInt ms_level) const;
    Iterator getClosestSpectrumInRT(const double RT, UInt ms_level);

    /// sets the chromatogram list
    void setChromatograms(const std::vector<MSChromatogram>& chromatograms);
    void setChromatograms(std::vector<MSChromatogram>&& chromatograms);

    /// adds a chromatogram to the list
    void addChromatogram(const MSChromatogram& chromatogram);
    void addChromatogram(MSChromatogram&& chrom);

    /// returns the chromatogram list
    const std::vector<MSChromatogram>& getChromatograms() const;

    /// returns the chromatogram list (mutable)
    std::vector<MSChromatogram>& getChromatograms();

    /// @name Easy Access interface
    //@{
    /// returns a single chromatogram
    MSChromatogram& getChromatogram(Size id);

    /// returns a single spectrum
    MSSpectrum& getSpectrum(Size id);

    /// get the total number of spectra available
    Size getNrSpectra() const;

    /// get the total number of chromatograms available
    Size getNrChromatograms() const;
    //@}

    /**
    @brief Computes the total ion chromatogram (TIC) for a given MS level (use ms_level = 0 for all levels). 

    By default, each MS spectrum's intensity just gets summed up. Regular RT bins can be obtained by specifying @p rt_bin_size.
    If a bin size in RT seconds greater than 0 is given resampling is used.

    @param rt_bin_size RT bin size in seconds (0 = no resampling)
    @param ms_level MS level of spectra for calculation (0 = all levels)
    @return TIC Chromatogram
    **/
    const MSChromatogram calculateTIC(float rt_bin_size = 0, UInt ms_level = 1) const;

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data);

    /// returns true if at least one of the spectra has the specified level
    bool containsScanOfLevel(size_t ms_level) const;

    /// returns true if any MS spectra of trthe specified level contain at least one peak with intensity of 0.0
    bool hasZeroIntensities(size_t ms_level) const;

    /// do any of the spectra have a PeptideID?
    bool hasPeptideIdentifications() const;

    /// Are all MSSpectra in this experiment part of an IM Frame? I.e. they all have the same RT, but different drift times
    bool isIMFrame() const;

  protected:
    /// MS levels of the data
    std::vector<UInt> ms_levels_;
    /// Number of all data points
    UInt64 total_size_;
    /// chromatograms
    std::vector<MSChromatogram > chromatograms_;
    /// spectra
    std::vector<SpectrumType> spectra_;

private:

    /// Helper class to add either general data points in set2DData or use mass traces from meta values
    template<typename ContainerValueType, bool addMassTraces>
    struct ContainerAdd_
    {
      static void addData_(SpectrumType* spectrum, const ContainerValueType* item);
      static void addData_(SpectrumType* spectrum, const ContainerValueType* item, const StringList& store_metadata_names);
    };

    template<typename ContainerValueType>
    struct ContainerAdd_<ContainerValueType, false>
    {
      /// general method for adding data points
      static void addData_(SpectrumType* spectrum, const ContainerValueType* item)
      {
        // create temporary peak and insert it into spectrum
        spectrum->insert(spectrum->end(), PeakType());
        spectrum->back().setIntensity(item->getIntensity());
        spectrum->back().setPosition(item->getMZ());
      }
      /// general method for adding data points, including metadata arrays (populated from metainfointerface)
      static void addData_(SpectrumType* spectrum, const ContainerValueType* item, const StringList& store_metadata_names)
      {
        addData_(spectrum, item);
        for (StringList::const_iterator itm = store_metadata_names.begin(); itm != store_metadata_names.end(); ++itm)
        {
          float val = std::numeric_limits<float>::quiet_NaN();
          if (item->metaValueExists(*itm)) val = item->getMetaValue(*itm);
          spectrum->getFloatDataArrays()[itm - store_metadata_names.begin()].push_back(val);
        }
      }
    };

    template<typename ContainerValueType>
    struct ContainerAdd_<ContainerValueType, true>
    {
      /// specialization for adding feature mass traces (does not support metadata_names currently)
      static void addData_(SpectrumType* spectrum, const ContainerValueType* item)
      {
        if (item->metaValueExists("num_of_masstraces"))
        {
          Size mts = item->getMetaValue("num_of_masstraces");
          int charge = (item->getCharge()==0 ? 1 : item->getCharge()); // set to 1 if charge is 0, otherwise div/0 below
          for (Size i = 0; i < mts; ++i)
          {
            String meta_name = String("masstrace_intensity_") + i;
            if (!item->metaValueExists(meta_name))
            {
              throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Meta value '") + meta_name + "' expected but not found in container.");
            }
            ContainerValueType p;
            p.setIntensity(item->getMetaValue(meta_name));
            p.setPosition(item->getMZ() + Constants::C13C12_MASSDIFF_U / charge * i);
            ContainerAdd_<ContainerValueType, false>::addData_(spectrum, &p);
          }
        }
        else ContainerAdd_<ContainerValueType, false>::addData_(spectrum, item);
      }
    };

    /*
      @brief Append a spectrum to current MSExperiment

      @param rt RT of new spectrum
      @return Pointer to newly created spectrum
    */
    SpectrumType* createSpec_(PeakType::CoordinateType rt);

    /*
      @brief Append a spectrum including floatdata arrays to current MSExperiment

      @param rt RT of new spectrum
      @param metadata_names Names of floatdata arrays attached to this spectrum
      @return Pointer to newly created spectrum
    */
    SpectrumType* createSpec_(PeakType::CoordinateType rt, const StringList& metadata_names);

  };

  /// Print the contents to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const MSExperiment& exp);

} // namespace OpenMS

#include <OpenMS/KERNEL/StandardTypes.h>



// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Samuel Wein $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>

#include <OpenMS/KERNEL/Peak1D.h>

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/DATASTRUCTURES/RegularExpression.h>

#include <algorithm>
#include <vector>

namespace OpenMS
{
  /**
      @brief Annotates spectra from identifications and theoretical spectra or
             identifications from spectra and theoretical spectra matching
             with various options

      @htmlinclude OpenMS_SpectrumAnnotator.parameters

      @ingroup Chemistry
  */
  class OPENMS_DLLAPI SpectrumAnnotator :
    public DefaultParamHandler
  {
    public:

    /** @name Constructors and Destructors
    */
    //@{
    /// default constructor
    SpectrumAnnotator();

    /// copy constructor
    SpectrumAnnotator(const SpectrumAnnotator & source);

    /// destructor
    ~SpectrumAnnotator() override;
    //@}

    /// assignment operator
    SpectrumAnnotator & operator=(const SpectrumAnnotator & tsg);

    /** @name Available annotators
     */
    //@{

    /**
      @brief Adds ion match annotation to the @p spec input spectrum

      @param[in] spec A PeakSpectrum containing the peaks from which the @p pi identifications are made
      @param[out] ph A spectrum identifications to be used for the annotation, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
      @param[in] tg A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
      @param[in] sa A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance

      The matches are added as DataArrays to the spectrum (Names: IonName and IonMatchError). The parameters of the TheoreticalSpectrumGenerator define the comprehensiveness of the available matching. The parameters of SpectrumAlignment define the matching tolerance.
      See the parameter section of this class for the available options.
    */
    void annotateMatches(PeakSpectrum &spec, const PeptideHit& ph, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const;

    /**
      @brief Adds ion match statistics to @p pi PeptideIdentifcation

      @param[in] pi A spectrum identifications to be annotated, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
      @param[in] spec A PeakSpectrum containing the peaks from which the @p pi identifications are made
      @param[in] tg A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
      @param[in] sa A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance

      The ion match statistics are added as UserParams to either the PeptideIdentification (parameters of the matching) and PeptideHit. The parameters of the TheoreticalSpectrumGenerator define the comprehensiveness of the available matching. The parameters of SpectrumAlignment define the matching tolerance.
      See the parameter section of this class for the available statistics.
    */
    void addIonMatchStatistics(PeptideIdentification& pi, MSSpectrum &spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const;

    /**
      @brief Adds peak annotations to the @p ph PeptideHit

      @param[in] ph A PeptideHit whose PeakAnnotations vector will be filled with the ion matches
      @param[in] spec A PeakSpectrum containing the peaks from which the @p ph identifications are made
      @param[in] tg A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
      @param[in] sa A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance
      @param[in] include_unmatched_peaks If true, all spectrum peaks will be included in the PeakAnnotations vector. 
             Unmatched peaks will have empty annotation strings. If false (default), only matched peaks are included.

      This method directly fills the PeakAnnotations vector of the PeptideHit with the matched ions.
      Each matched peak gets a PeakAnnotation containing the annotation string (ion name), charge,
      observed m/z, intensity, and theoretical m/z.
      The parameters of the TheoreticalSpectrumGenerator define the comprehensiveness of the available matching.
      The parameters of SpectrumAlignment define the matching tolerance.
    */
    void addPeakAnnotationsToPeptideHit(PeptideHit& ph, const PeakSpectrum& spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa, bool include_unmatched_peaks = false) const;

    /**
      @brief Count matched rungs in the longest connected ion series.

      Positions refer to one ion series, independently of molecule type, charge,
      and neutral-loss annotation. A missing rung does not contribute to the
      returned count. With the default gap limit, only consecutive matches are
      connected. Callers must keep distinct ion series and fragment models
      separate before calling this function.

      @param[in] matched_positions Whether each ordinal position was matched.
      @param[in] max_missing Maximum consecutive missing positions allowed between matches.
      @return Number of matched positions in the longest connected series, or zero.
    */
    static Size computeLongestIonSeries(const std::vector<bool>& matched_positions, Size max_missing = 0)
    {
      Size longest = 0, current = 0, missing = 0;
      for (bool matched : matched_positions)
      {
        if (matched)
        {
          ++current;
          missing = 0;
          longest = std::max(longest, current);
        }
        else if (current != 0 && ++missing > max_missing) { current = 0; }
      }
      return longest;
    }

    /**
      @brief Count connected matched ordinals without materializing a position mask.

      Ordinals must be sorted in ascending order. Duplicate positions contribute
      once, even when several charges or peaks match the same rung. The caller
      supplies a projection so existing annotation or match records can be used
      directly without copying them.

      @param[in] first First match record.
      @param[in] last Past-the-end match record.
      @param[in] max_missing Maximum consecutive missing positions between matches.
      @param[in] position Projection from a match record to its ordinal position.
      @return Number of distinct matched ordinals in the longest connected series.
    */
    template<typename Iterator, typename Position>
    static Size computeLongestIonSeries(Iterator first, Iterator last, Size max_missing, Position position)
    {
      Size longest = 0, current = 0, previous = 0;
      for (; first != last; ++first)
      {
        const Size ordinal = position(*first);
        OPENMS_PRECONDITION(current == 0 || ordinal >= previous, "Ion ordinals must be sorted");
        if (current != 0 && ordinal == previous) { continue; }
        current = current != 0 && ordinal - previous - 1 <= max_missing ? current + 1 : 1;
        previous = ordinal;
        longest = std::max(longest, current);
      }
      return longest;
    }

    /// overwrite
    void updateMembers_() override;

    //@}

    protected:
      bool basic_statistics_;
      bool list_of_ions_matched_;
      bool max_series_;
      bool SN_statistics_;
      bool precursor_statistics_;
      unsigned topNmatch_fragmenterrors_;
      bool fragmenterror_statistics_;
      bool terminal_series_match_ratio_;

      static const RegularExpression nt_regex_;
      static const RegularExpression ct_regex_;
      static const RegularExpression noloss_regex_;
      static const RegularExpression seriesposition_regex_;
  };
}

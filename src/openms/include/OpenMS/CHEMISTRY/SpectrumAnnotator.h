// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>

#include <OpenMS/KERNEL/Peak1D.h>

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <boost/regex.hpp>

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

      @param spec A PeakSpectrum containing the peaks from which the @p pi identifications are made
      @param ph A spectrum identifications to be used for the annotation, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
      @param tg A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
      @param sa A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance

      The matches are added as DataArrays to the spectrum (Names: IonName and IonMatchError). The parameters of the TheoreticalSpectrumGenerator define the comprehensiveness of the available matching. The parameters of SpectrumAlignment define the matching tolerance.
      See the parameter section of this class for the available options.
    */
    void annotateMatches(PeakSpectrum &spec, const PeptideHit& ph, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const;

    /**
      @brief Adds ion match statistics to @p pi PeptideIdentifcation

      @param pi A spectrum identifications to be annotated, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
      @param spec A PeakSpectrum containing the peaks from which the @p pi identifications are made
      @param tg A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
      @param sa A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance

      The ion match statistics are added as UserParams to either the PeptideIdentification (parameters of the matching) and PeptideHit. The parameters of the TheoreticalSpectrumGenerator define the comprehensiveness of the available matching. The parameters of SpectrumAlignment define the matching tolerance.
      See the parameter section of this class for the available statistics.
    */
    void addIonMatchStatistics(PeptideIdentification& pi, MSSpectrum &spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const;

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

      static const boost::regex nt_regex_;
      static const boost::regex ct_regex_;
      static const boost::regex noloss_regex_;
      static const boost::regex seriesposition_regex_;
  };
}


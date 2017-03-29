// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_SPECTRUMANNOTATOR_H
#define OPENMS_CHEMISTRY_SPECTRUMANNOTATOR_H

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <boost/regex.hpp>

namespace OpenMS
{
  /**
      @brief Generates theoretical spectra with various options

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
    virtual ~SpectrumAnnotator();
    //@}

    /// assignment operator
    SpectrumAnnotator & operator=(const SpectrumAnnotator & tsg);

    /** @name Available annotators
     */
    //@{

    /**
      @brief Generates an ion match annotated version of the @p spec input spectrum

      @param ph A spectrums identifications to be used for the annotation, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
      @param spec A PeakSpectrum containing the peaks from which the @p pi identifications are made
      @param tg A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
      @param sa A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance

      The matches are added as UserParams to the peaks (MetaValues IonName and IonMatchError). The parameters of the TheoreticalSpectrumGenerator define the comprehensiveness of the available matching. The parameters of SpectrumAlignment define the matching tolerance.
      See the parameter section of this class for the available options.
      @htmlinclude OpenMS_SpectrumAnnotator.parameters
    */
    MSSpectrum<RichPeak1D> annotateMatches(const PeptideHit& ph, const MSSpectrum<Peak1D>& spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const;

    /**
      @brief Adds ion match statistics to @p pi PeptideIdentifcation

      @param pi A spectrums identifications to be annotated, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
      @param spec A PeakSpectrum containing the peaks from which the @p pi identifications are made
      @param tg A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
      @param sa A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance

      The ion match statistics are added as UserParams to either the PeptideIdentification (parameters of the matching) and PeptideHit. The parameters of the TheoreticalSpectrumGenerator define the comprehensiveness of the available matching. The parameters of SpectrumAlignment define the matching tolerance.
      See the parameter section of this class for the available statistics.
      @htmlinclude OpenMS_SpectrumAnnotator.parameters
    */
    void addIonMatchStatistics(PeptideIdentification& pi, const MSSpectrum<Peak1D>& spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const;

    /// overwrite
    void updateMembers_();

    //@}

    protected:
      bool basic_statistics_;
      bool list_of_ions_matched_;
      bool max_series_;
      bool SN_statistics_;
      bool precursor_statistics_;
      uint topNmatch_fragmenterrors_;
      bool fragmenterror_statistics_;
      bool terminal_series_match_ratio_;

      static const boost::regex nt_regex_;
      static const boost::regex ct_regex_;
      static const boost::regex noloss_regex_;

  };
}
#endif


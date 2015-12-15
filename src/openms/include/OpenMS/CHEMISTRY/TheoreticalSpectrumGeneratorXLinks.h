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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef THEORETICALSPECTRUMGENERATORXLINKS_H
#define THEORETICALSPECTRUMGENERATORXLINKS_H

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>


namespace OpenMS
{
  class AASequence;

  /**
      @brief Generates theoretical spectra with various options

  @htmlinclude OpenMS_TheoreticalSpectrumGenerator.parameters

      @ingroup Chemistry
  */
  class OPENMS_DLLAPI TheoreticalSpectrumGeneratorXLinks :
    public DefaultParamHandler
  {
public:

    struct ProteinProteinCrossLink
    {
      /** Enums
      */
      //@{
      /** @brief type of Protein-Protein cross-link
      */
      enum ProteinProteinCrossLinkType
      {
        PROTEIN_PROTEIN = 0,
        MONO = 1,
        LOOP = 2,
        NUMBER_OF_TERM_SPECIFICITY
      };

      AASequence alpha; // longer peptide
      AASequence beta; // shorter peptide (empty for mono-link), tie bracker: mass then lexicographical
      std::pair<Size, Size> cross_link_position; // index in alpha, beta or between alpha, alpha in mono-links

      ProteinProteinCrossLinkType getType()
      {
        if (!beta.empty()) return PROTEIN_PROTEIN;

        if (cross_link_position.second == 0) return MONO;

        return LOOP;
      }

      double getMass(double cross_linker_mass)
      {
        switch(getType())
        {
          case PROTEIN_PROTEIN: return 0; break;
          case MONO: return 0; break;
          case LOOP: return 0; break;
          default: 
          //TODO: error
          break;
        }
        return cross_linker_mass; // TODO change
      }
    };


    /** @name Constructors and Destructors
    */
    //@{
    /// default constructor
    TheoreticalSpectrumGeneratorXLinks();

    /// copy constructor
    TheoreticalSpectrumGeneratorXLinks(const TheoreticalSpectrumGeneratorXLinks & source);

    /// destructor
    virtual ~TheoreticalSpectrumGeneratorXLinks();
    //@}

    /// assignment operator
    TheoreticalSpectrumGeneratorXLinks & operator=(const TheoreticalSpectrumGeneratorXLinks & tsg);

    /** @name Acessors
     */
    //@{
    /// returns a spectrum with b and y peaks
    virtual void getSpectrum(RichPeakSpectrum & spec, const ProteinProteinCrossLink & cross_link, Int charge = 1) const;
    virtual void getCommonIonSpectrum(RichPeakSpectrum & spec, const ProteinProteinCrossLink & cross_link, Int charge = 1) const;
    virtual void getCommonIonSpectrum(RichPeakSpectrum & spec, const AASequence & peptide, Int charge = 1) const;
    virtual void getXLinkIonSpectrum(RichPeakSpectrum & spec_alpha, RichPeakSpectrum & spec_beta, const ProteinProteinCrossLink& cross_link, Int mincharge = 3, Int maxcharge = 5) const;

    /// adds peaks to a spectrum of the given ion-type, peptide, charge, and intensity
    virtual void addPeaks(RichPeakSpectrum & spectrum, const ProteinProteinCrossLink & cross_link, Residue::ResidueType res_type, Int charge = 1) const;
    virtual void addCommonPeaks(RichPeakSpectrum & spectrum, const AASequence & peptide, Residue::ResidueType res_type, Int charge = 1) const;
    virtual void addXLinkIonPeaks(RichPeakSpectrum & spec_alpha, RichPeakSpectrum & spec_beta, const ProteinProteinCrossLink & cross_link, Residue::ResidueType res_type, Int charge) const;

    /// adds the precursor peaks to the spectrum
    virtual void addPrecursorPeaks(RichPeakSpectrum & spec, const AASequence & peptide, Int charge = 1) const;

    /// Adds the common, most abundant immonium ions to the theoretical spectra if the residue is contained in the peptide sequence
    void addAbundantImmoniumIons(RichPeakSpectrum & spec, const AASequence& peptide) const;
    //@}
  };
}

#endif // THEORETICALSPECTRUMGENERATORXLINKS_H

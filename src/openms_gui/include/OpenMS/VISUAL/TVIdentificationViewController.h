// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <OpenMS/VISUAL/TVControllerBase.h>
#include <vector>

namespace OpenMS
{
  class NASequence;
  class SpectraIDViewTab;
  class TOPPViewBase;

  /**
  @brief Behavior of TOPPView in identification mode.
  */
  class TVIdentificationViewController
    : public TVControllerBase
  {
    Q_OBJECT

  public:
    /// Construct the behaviour with its parent
    TVIdentificationViewController(TOPPViewBase* parent, SpectraIDViewTab* spec_id_view_);

  public slots:
    /// Behavior for showSpectrumAsNew1D
    virtual void showSpectrumAsNew1D(int spectrum_index, int peptide_id_index, int peptide_hit_index);

    /// Show spectrum without selecting an identification
    virtual void showSpectrumAsNew1D(int index);

    /// Behavior for activate1DSpectrum
    virtual void activate1DSpectrum(int spectrum_index, int peptide_id_index, int peptide_hit_index);

    /// select spectrum without selecting an identification
    virtual void activate1DSpectrum(int index);

    /// Behavior for deactivate1DSpectrum
    virtual void deactivate1DSpectrum(int index);

    /// Slot for behavior activation
    void activateBehavior() override;

    /// Slot for behavior deactivation
    void deactivateBehavior() override;

    void setVisibleArea1D(double l, double h);

  private:
    /// Adds labels for the provided precursors to the 1D spectrum
    void addPrecursorLabels1D_(const std::vector<Precursor>& pcs);

    /// Removes the precursor labels for from the specified 1D spectrum
    void removeTemporaryAnnotations_(Size spectrum_index);

    /// Adds a theoretical spectrum as set from the preferences dialog for the peptide hit.
    void addTheoreticalSpectrumLayer_(const PeptideHit& ph);

    /// Add peak annotatios from id data structure
    void addPeakAnnotationsFromID_(const PeptideHit& hit);

    /// removes all layer with theoretical spectrum generated in identification view
    void removeTheoreticalSpectrumLayer_();

    /// remove all graphical peak annotations
    void removeGraphicalPeakAnnotations_(int spectrum_index);

    /// Adds annotation (compound name, adducts, ppm error) to a peak in 1D spectra
    void addPeakAnnotations_(const std::vector<PeptideIdentification>& ph);

    /// Helper function for text formatting
    String n_times(Size n, String input);

    /// Helper function that turns fragment annotations into coverage Strings for visualization with the sequence
    void extractCoverageStrings(std::vector<PeptideHit::PeakAnnotation> frag_annotations, String& alpha_string, String& beta_string, Size alpha_size, Size beta_size);

    /// Generates HTML for showing the sequence with annotations of matched fragments
    template <typename SeqType>
    String generateSequenceDiagram_(const SeqType& seq, const std::vector<PeptideHit::PeakAnnotation>& annotations, const StringList& top_ions, const StringList& bottom_ions);

    /// Helper function for generateSequenceDiagram_() - overload for peptides
    void generateSequenceRow_(const AASequence& seq, std::vector<String>& row);

    /// Helper function for generateSequenceDiagram_() - overload for oligonucleotides
    void generateSequenceRow_(const NASequence& seq, std::vector<String>& row);

    /// Helper function, that collapses a vector of Strings into one String
    String collapseStringVector(std::vector<String> strings);

  private:
    SpectraIDViewTab* spec_id_view_;
    /// Used to check which annotation handles have been added automatically by the identification view. Ownership
    /// of the AnnotationItems has the Annotation1DContainer
    std::vector<Annotation1DItem*> temporary_annotations_;
  };
}

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

#include <QObject>

#include <OpenMS/METADATA/SpectrumSettings.h>

#include <vector>

namespace OpenMS
{
  class TopDownViewWidget;
  class TOPPViewBase;
  class Annotation1DItem;

  /**
  @brief Behavior of TOPPView in identification mode.
  */
  class TOPPViewTopDownViewBehavior
    : public QObject
  {
    Q_OBJECT

  public:
    /// Construct the behaviour with its parent
    TOPPViewTopDownViewBehavior(TOPPViewBase* parent, TopDownViewWidget* spec_id_view_);

  public slots:
    /// Show spectrum without selecting an identification
    virtual void showSpectrumAs1D(int index);

    /// select spectrum without selecting an identification
    virtual void activate1DSpectrum(int index);

    /// Behavior for deactivate1DSpectrum
    virtual void deactivate1DSpectrum(int index);

    /// Slot for behavior activation
    virtual void activateBehavior();

    /// Slot for behavior deactivation
    virtual void deactivateBehavior();

    /// Slot to change visible area to [l, h] in m/z / mass coordinates
    void setVisibleArea1D(double l, double h);

  private:
    /// Adds labels for the provided precursors to the 1D spectrum
    void addPrecursorLabels1D_(const std::vector<Precursor>& pcs);

    /// Removes the precursor labels for from the specified 1D spectrum
    void removeTemporaryAnnotations_(Size spectrum_index);

    /// remove all graphical peak annotations
    void removeGraphicalPeakAnnotations_(int spectrum_index);

    /// annotate fixed masses TODO: add code
    void addPeakAnnotations_(const std::vector<double>& masses);
  private:
    TOPPViewBase* tv_;
    TopDownViewWidget* spec_id_view_;
    /// Used to check which annotation handles have been added automatically by the identification view. Ownership
    /// of the AnnotationItems has the Annotation1DContainer
    std::vector<Annotation1DItem*> temporary_annotations_;
  };
}

// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Tom Waschischeck $
// --------------------------------------------------------------------------


#pragma once
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtWidgets/QDialog>
#include <QtWidgets/qspinbox.h>
#include <QtWidgets/qlabel.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <array>

namespace Ui
{
  class TheoreticalSpectrumGenerationDialogTemplate;
}

class QListWidgetItem;

namespace OpenMS
{
  class TestTSGDialog; // fwd declaring test class

  /// state of an ion (and its intensity)
  enum class CheckBoxState
  {
    HIDDEN,    ///< check box hidden (invisible)
    ENABLED,   ///< check box enabled (visible, but not checked)
    PRECHECKED ///< check box enabled and checked by default
  };

  /**
      @brief Dialog which allows to enter an AA or NA sequence and generates a theoretical spectrum for it.

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TheoreticalSpectrumGenerationDialog :
    public QDialog
  {
    Q_OBJECT

  public:
    /// struct for all information about a check box of an ion
    struct CheckBox {
      /// Constructor
      CheckBox(QDoubleSpinBox** sb, QLabel** l, std::array<CheckBoxState, 3> s, std::pair<String, String> p_t, std::pair<String, String> p_s);

      /// pointer to the corresponding ion intensity spin box
      QDoubleSpinBox** ptr_to_spin_box;

      /// pointer to the label of the spin box
      QLabel** ptr_to_spin_label;

      /// State of this check box depending on sequence type ("Peptide", "RNA", "Metabolite")
      const std::array<CheckBoxState, 3> state;

      /// parameter with description of this ion
      const std::pair<String, String> param_this;

      /// parameter with description of the ion intensity
      const std::pair<String, String> param_spin;
    };

    /// type of the input sequence (corresponds to the value of the combo box 'ui_->seq_type')
    enum class SequenceType
    {
      PEPTIDE,
      RNA,
      METABOLITE
    };
    
    friend class TestTSGDialog; // to test the GUI expressed in the private member ui

    /// Constructor
    TheoreticalSpectrumGenerationDialog();
    /// Destructor
    ~TheoreticalSpectrumGenerationDialog() override;

    /// returns the calculated spectrum
    const MSSpectrum& getSpectrum() const;

    /// returns the input sequence (is public for TOPPView)
    const String getSequence() const;

protected slots:

    /// for isotope model changes
    void modelChanged_();
    /// for sequence type changes (combo box)
    void seqTypeSwitch_();
    /// change check state of check box on widget click
    void listWidgetItemClicked_(QListWidgetItem* item);
    /// calculates the spectrum
    void calculateSpectrum_();

protected:

private:
    /// calculate parameters from UI elements
    Param getParam_() const;

    /// iterates through 'check_boxes_' and en-/disables
    /// check boxes and corresponding spin boxes (and their labels)
    void updateIonTypes_();

    /// UI
    Ui::TheoreticalSpectrumGenerationDialogTemplate* ui_;

    /// save current sequence setting
    SequenceType seq_type_;

    /// array of TSGDialog::CheckBox
    /// 
    /// Note: Ordering has to be the same as in the UI!
    const std::array<CheckBox, 12> check_boxes_;

    /// member to save the calculated spectrum to
    MSSpectrum spec_;
  };

}

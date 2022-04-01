// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

  // Note: If an additional check box is added all of the following three objects have to be edited!
  //
  //
  // enum to get ion check box index (Ordering has to be the same as in the ui!)
  enum class CheckBox
  {
    A_Ions,
    A_b_Ions,
    B_Ions,
    C_Ions,
    D_Ions,
    W_Ions,
    X_Ions,
    Y_Ions,
    Z_Ions,
    Precursor,
    Neutral_losses,
    Abundant_Immonium_Ions,
    NUMBER_OF_CHECK_BOXES
  };

  enum class State
  {
    HIDDEN,
    ENABLED,
    PRECHECKED
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
    struct CheckBox_ {
      CheckBox_(QDoubleSpinBox** p, std::array<State, 3> s, std::pair<String, String> p_t, std::pair<String, String> p_s);

      QDoubleSpinBox** ptr_to_spin_box;

      std::array<State, 3> state;

      std::pair<String, String> param_this;

      std::pair<String, String> param_spin;
    };
    
    friend class TestTSGDialog; // to test the GUI expressed in the private member ui

    /// Constructor
    TheoreticalSpectrumGenerationDialog();
    /// Destructor
    ~TheoreticalSpectrumGenerationDialog() override;

    const MSSpectrum& getSpectrum() const;

    const String getSequence() const;

protected slots:

    void modelChanged_();
    void calculateSpectrum_();
    void seqTypeSwitch_();
    void listWidgetItemClicked_(QListWidgetItem* item);

protected:

private:
    QString seq_type_;

    Param getParam_() const;

    Ui::TheoreticalSpectrumGenerationDialogTemplate* ui_;

    MSSpectrum spec_;

    const std::array<CheckBox_, 12> check_boxes_;

    // For each check box get its intensity spin box and its ion name aka parameter prefix.
    // Order is important here!
    // To access the right entry for each check box
    // use int(TheoreticalSpectrumGenerationDialog::CheckBox).
    const std::vector<std::pair<QDoubleSpinBox**, String>> check_box_to_intensity_;

    // To check if the check box is enabled for <'Peptide','RNA'>
    // Note that for metabolite input there are no ions.
    // Order is important here!
    // To access the right entry for each check box
    // use int(TheoreticalSpectrumGenerationDialog::CheckBox).
    const std::vector<std::pair<bool, bool>> intensity_ion_exists {
      {1, 1}, // A-Ion
      {0, 1}, // a-B-Ion
      {1, 1}, // B-Ion
      {1, 1}, // C-Ion
      {0, 1}, // D-Ion
      {0, 1}, // W-Ion
      {1, 1}, // X-Ion
      {1, 1}, // Y-Ion
      {1, 1}, // Z-Ion
      {1, 1}, // Precursor
      {1, 0}, // Neutral losses
      {1, 0}  // Abundant Immonium Ions
    };
  };

}

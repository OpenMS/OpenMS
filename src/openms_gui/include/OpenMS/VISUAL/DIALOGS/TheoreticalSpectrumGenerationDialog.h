// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

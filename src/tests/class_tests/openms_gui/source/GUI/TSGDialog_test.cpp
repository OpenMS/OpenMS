// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <TSGDialog_test.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <QtWidgets/QPushButton>
#include <QtWidgets/QMessageBox>

#include <qdebug.h>

using namespace OpenMS;
using namespace std;

// delay in ms
// higher values together with 'dialog_.show();' can be useful for debugging this test
constexpr int DELAY {5};

template<typename T>
void TestTSGDialog::testSpinBox_(T* box, string str_value)
{
  if (box == nullptr) return;
  
  double value = stod(str_value);

  // skip illegal input
  if (!(value >= box->minimum() && value <= box->maximum())) return;
  
  box->clear();

  // simulate keyboard input
  QTest::keyClicks(box, QString::fromStdString(str_value));
  QTest::qWait(DELAY);

  // verify change
  QVERIFY(value == double(box->value())); // double cast needed because of template
}

void TestTSGDialog::testIsotopeModel_(bool skip_none)
{
  // QTest::mouseClick needs the exact position of the interactable part of the radio button
  // If someone knows how to get this position, please adapt this code.

  if (skip_none)
  {
    QVERIFY(!UI->model_none_button->isEnabled());
  }
  else
  {
    UI->model_none_button->click();
    QTest::qWait(DELAY);
    QVERIFY(!(UI->max_iso_spinbox->isEnabled()));
    QVERIFY(!(UI->max_iso_label->isEnabled()));
    QVERIFY(!(UI->max_iso_prob_spinbox->isEnabled()));
    QVERIFY(!(UI->max_iso_prob_label->isEnabled()));
  }

  UI->model_coarse_button->click();
  QTest::qWait(DELAY);
  QVERIFY(UI->max_iso_spinbox->isEnabled());
  QVERIFY(UI->max_iso_label->isEnabled());
  QVERIFY(!(UI->max_iso_prob_spinbox->isEnabled()));
  QVERIFY(!(UI->max_iso_prob_label->isEnabled()));
  testSpinBox_(UI->max_iso_spinbox);

  UI->model_fine_button->click();
  QTest::qWait(DELAY);
  QVERIFY(!(UI->max_iso_spinbox->isEnabled()));
  QVERIFY(!(UI->max_iso_label->isEnabled()));
  QVERIFY(UI->max_iso_prob_spinbox->isEnabled());
  QVERIFY(UI->max_iso_prob_label->isEnabled());
  testSpinBox_(UI->max_iso_prob_spinbox);
}

void OpenMS::TestTSGDialog::testIonsIntensities_()
{
  int input_type;
  if (dialog_.seq_type_ == TheoreticalSpectrumGenerationDialog::SequenceType::PEPTIDE)
    input_type = 0;
  else if (dialog_.seq_type_ == TheoreticalSpectrumGenerationDialog::SequenceType::RNA)
    input_type = 1;
  else // Metabolite
    input_type = 2;

  for (size_t i = 0; i < dialog_.check_boxes_.size(); ++i)
  {
    // get the item
    QListWidgetItem* item = UI->ion_types->item(i);
    QVERIFY(item);

    const TheoreticalSpectrumGenerationDialog::CheckBox* curr_box = &dialog_.check_boxes_.at(i);

    // get intensity spin box corresponding to current check box
    QDoubleSpinBox** spin_ptr = curr_box->ptr_to_spin_box;

    if (curr_box->state.at(input_type) != CheckBoxState::HIDDEN)
    {
      // check state before clicking
      Qt::CheckState prev = item->checkState();

      // get the rectangular coordinates of the item
      QRect rect = UI->ion_types->visualItemRect(item);

      // imitate the click on check box c
      QTest::mouseClick(UI->ion_types->viewport(), Qt::LeftButton, Qt::NoModifier, rect.center());
      QTest::qWait(DELAY);

      // verfiy the check state changed
      QVERIFY2(prev != item->checkState(), qPrintable("Clicking on '" + item->text() + "' didn't change its state."));

      if (spin_ptr == nullptr)
        continue;
      testSpinBox_(*spin_ptr);
    }
    else
    {
      // if ion type isn't supported, check if the ion and its intensity are hidden
      QVERIFY2(item->isHidden(), qPrintable(item->text() + " wasn't hidden."));
      if (spin_ptr == nullptr)
        continue;
      QVERIFY2((*spin_ptr)->isHidden(), qPrintable("Spin box of '" + item->text() + "' wasn't hidden."));
    }
  }
}

void OpenMS::TestTSGDialog::testSequenceInput_(QString input)
{
  UI->seq_input->clear();
  QTest::keyClicks(UI->seq_input, input);
  QString read_seq = QString::fromStdString(std::string(dialog_.getSequence()));
  QVERIFY2(read_seq == UI->seq_input->text(), "Test of sequence input failed!");
}

void OpenMS::TestTSGDialog::checkMessageBoxExists_()
{
  // get the active window
  QWidget* active_widget = QApplication::activeModalWidget();
  if (active_widget->inherits("QMessageBox")) // if it's a message box, close it
  {
    QMessageBox* mb = qobject_cast<QMessageBox*>(active_widget);
    QTest::keyClick(mb, Qt::Key_Enter);
    QVERIFY(true);
    return;
  }
  QVERIFY2(false, "Expected message box wasn't produced!");
}

void TestTSGDialog::testMessageBoxes_()
{
  // check empty sequence input
  UI->seq_input->clear();
  QTimer::singleShot(DELAY, this, &TestTSGDialog::checkMessageBoxExists_); // calls the SLOT after DELAY ms passed
  UI->dialog_buttons->button(QDialogButtonBox::Ok)->click();

  // check unparsible sequence input
  UI->seq_input->setText("J+-5!");
  QTimer::singleShot(DELAY, this, &TestTSGDialog::checkMessageBoxExists_);
  UI->dialog_buttons->button(QDialogButtonBox::Ok)->click();

  // unselect all ions to produce empty spectrum
  for (size_t i = 0; i < dialog_.check_boxes_.size(); ++i)
  {
    UI->ion_types->item(i)->setCheckState(Qt::CheckState::Unchecked);
  }
  QTimer::singleShot(DELAY, this, &TestTSGDialog::checkMessageBoxExists_);
  UI->dialog_buttons->button(QDialogButtonBox::Ok)->click();
}

void TestTSGDialog::testConstruction()
{
  // editable/interactable GUI parts
  QVERIFY2(UI->seq_type, "Sequence selection combo box not created.");
  QVERIFY2(UI->seq_input, "Sequence input line edit not created.");
  QVERIFY2(UI->charge_spinbox, "Charge spin box not created.");
  QVERIFY2(UI->max_iso_spinbox, "Max. isotope model spin box not created.");
  QVERIFY2(UI->max_iso_prob_spinbox, "Max. isotope probability spin box not created.");
  QVERIFY2(UI->ion_types, "Ion list widget not created.");
  QVERIFY2(UI->a_intensity, "A ion intensity spin box not created.");
  QVERIFY2(UI->a_b_intensity, "A-b ion intensity spin box not created.");
  QVERIFY2(UI->b_intensity, "B ion intensity spin box not created.");
  QVERIFY2(UI->c_intensity, "C ion intensity spin box not created.");
  QVERIFY2(UI->d_intensity, "D ion intensity spin box not created.");
  QVERIFY2(UI->w_intensity, "W ion intensity spin box not created.");
  QVERIFY2(UI->x_intensity, "X ion intensity spin box not created.");
  QVERIFY2(UI->y_intensity, "Y ion intensity spin box not created.");
  QVERIFY2(UI->z_intensity, "Z ion intensity spin box not created.");
  QVERIFY2(UI->rel_loss_intensity, "Relative loss intensity spin box not created.");
  QVERIFY2(UI->dialog_buttons, "Buttonbox not created.");

  // labels
  QVERIFY2(UI->enter_seq_label, "'Enter sequence' label not created.");
  QVERIFY2(UI->charge_label, "'Charge' label not created.");
  QVERIFY2(UI->ion_types_label, "'Generate' label not created.");
  QVERIFY2(UI->max_iso_label, "'Max. Isotope' label not created.");
  QVERIFY2(UI->max_iso_prob_label, "'Max. Isotope Probability in %' label not created.");
  QVERIFY2(UI->a_label, "'A-ions' label not created.");
  QVERIFY2(UI->a_b_label, "'A-b-ions' label not created.");
  QVERIFY2(UI->b_label, "'B-ions' label not created.");
  QVERIFY2(UI->c_label, "'C-ions' label not created.");
  QVERIFY2(UI->d_label, "'D-ions' label not created.");
  QVERIFY2(UI->w_label, "'W-ions' label not created.");
  QVERIFY2(UI->x_label, "'X-ions' label not created.");
  QVERIFY2(UI->y_label, "'Y-ions' label not created.");
  QVERIFY2(UI->z_label, "'Z-ions' label not created.");
  QVERIFY2(UI->rel_loss_label, "'Relative loss in %' label not created.");
  
  // group boxes
  QVERIFY2(UI->isotope_model, "Isotope model group box not created.");
  QVERIFY2(UI->intensities, "Intensity group box not created.");
}

void TestTSGDialog::testGui()
{
  //////////////////////////////////////////////////////
  //                     PEPTIDE                      //
  //////////////////////////////////////////////////////
  UI->seq_type->setCurrentText("Peptide");
  QTest::qWait(DELAY);

  // sequence input
  testSequenceInput_("PEPTIDE");

  // charge
  testSpinBox_(UI->charge_spinbox);

  // isotope model
  QVERIFY2(!UI->isotope_model->isHidden(), "Isotope model hidden for 'Peptide' setting.");
  testIsotopeModel_();

  // ion types and intensities
  testIonsIntensities_();

  // check relative loss intensity manually
  UI->rel_loss_intensity->clear();
  QTest::keyClicks(UI->rel_loss_intensity, "2");
  QTest::qWait(DELAY);
  QVERIFY(2.0 == UI->rel_loss_intensity->value());

  //////////////////////////////////////////////////////
  //                      RNA                         //
  //////////////////////////////////////////////////////
  UI->seq_type->setCurrentText("RNA");
  QTest::qWait(DELAY);

  // sequence input
  testSequenceInput_("ACGUGCA");

  // charge
  testSpinBox_(UI->charge_spinbox);

  // isotope model
  QVERIFY2(UI->isotope_model->isHidden(), "Isotope model not hidden for 'Peptide' setting.");

  // ion types and intensities
  testIonsIntensities_();
  // check relative loss intensity manually
  QVERIFY(UI->rel_loss_intensity->isHidden());
  QVERIFY(UI->rel_loss_label->isHidden());

  //////////////////////////////////////////////////////
  //                   Metabolite                     //
  //////////////////////////////////////////////////////
  UI->seq_type->setCurrentText("Metabolite");
  QTest::qWait(DELAY);

  // sequence input
  testSequenceInput_("C100H70N2O6");

  // charge
  testSpinBox_(UI->charge_spinbox);

  // isotope model
  QVERIFY2(!UI->isotope_model->isHidden(), "Isotope model hidden for 'Metabolite' setting.");
  testIsotopeModel_(true);

  // ion types and intensities
  //testIonsIntensities_();
}

void TestTSGDialog::testParameterImport()
{
  // set some parameters
  UI->seq_type->setCurrentText("Peptide");
  UI->charge_spinbox->setValue(3);
  UI->model_coarse_button->click();
  UI->max_iso_spinbox->setValue(5);
  for (size_t i = 0; i < dialog_.check_boxes_.size(); ++i)
  {
    UI->ion_types->item(i)->setCheckState(Qt::CheckState::Checked); // just check all boxes

    // get intensity spin box corresponding to current check box
    QDoubleSpinBox** spin_ptr = dialog_.check_boxes_.at(i).ptr_to_spin_box;
    if (spin_ptr == nullptr) continue;

    (*spin_ptr)->setValue(1.23);
  }
  UI->rel_loss_intensity->setValue(16);

  // get the parameters from the dialog
  Param p = dialog_.getParam_();

  // check if the returned parameters are correct
  QVERIFY2(int(p.getValue("charge")) == 3, "Parameter 'charge' wasn't set correctly.");
  QVERIFY2(string(p.getValue("isotope_model")) == "coarse", "Parameter 'isotope_model' wasn't set correctly. Expected 'coarse'.");
  QVERIFY2(int(p.getValue("max_isotope")) == 5, "Parameter 'max_isotope' wasn't set correctly.");
  QVERIFY2(string(p.getValue("add_a_ions")) == "true", "Parameter 'add_a_ions' wasn't set correctly.");
  QVERIFY2(double(p.getValue("a_intensity")) == 1.23, "Parameter 'a_intensity' wasn't set correctly.");

  QVERIFY2(string(p.getValue("add_b_ions")) == "true", "Parameter 'add_b_ions' wasn't set correctly.");
  QVERIFY2(double(p.getValue("b_intensity")) == 1.23, "Parameter 'b_intensity' wasn't set correctly.");

  QVERIFY2(string(p.getValue("add_c_ions")) == "true", "Parameter 'add_c_ions' wasn't set correctly.");
  QVERIFY2(double(p.getValue("c_intensity")) == 1.23, "Parameter 'c_intensity' wasn't set correctly.");

  QVERIFY2(string(p.getValue("add_x_ions")) == "true", "Parameter 'add_x_ions' wasn't set correctly.");
  QVERIFY2(double(p.getValue("x_intensity")) == 1.23, "Parameter 'x_intensity' wasn't set correctly.");

  QVERIFY2(string(p.getValue("add_y_ions")) == "true", "Parameter 'add_y_ions' wasn't set correctly.");
  QVERIFY2(double(p.getValue("y_intensity")) == 1.23, "Parameter 'y_intensity' wasn't set correctly.");

  QVERIFY2(string(p.getValue("add_z_ions")) == "true", "Parameter 'add_z_ions' wasn't set correctly.");
  QVERIFY2(double(p.getValue("z_intensity")) == 1.23, "Parameter 'z_intensity' wasn't set correctly.");

  QVERIFY2(string(p.getValue("add_losses")) == "true", "Parameter 'add_losses' wasn't set correctly.");
  QVERIFY2(double(p.getValue("relative_loss_intensity")) == 0.16, "Parameter 'relative_loss_intensity' wasn't set correctly.");

  QVERIFY2(string(p.getValue("add_abundant_immonium_ions")) == "true", "Parameter 'add_abundant_immonium_ions' wasn't set correctly.");
  QVERIFY2(string(p.getValue("add_precursor_peaks")) == "true", "Parameter 'add_precursor_peaks' wasn't set correctly.");

  // try the other isotope models
  UI->model_none_button->click();
  p.clear();
  p = dialog_.getParam_();
  QVERIFY2(string(p.getValue("isotope_model")) == "none", "Parameter 'isotope_model' wasn't set correctly. Expected 'none'.");

  UI->model_fine_button->click();
  UI->max_iso_prob_spinbox->setValue(16);
  p.clear();
  p = dialog_.getParam_();
  QVERIFY2(string(p.getValue("isotope_model")) == "fine", "Parameter 'isotope_model' wasn't set correctly. Expected 'fine'.");
  QVERIFY2(double(p.getValue("max_isotope_probability")) == 0.16, "Parameter 'max_isotope_probability' wasn't set correctly.");

  // check remaining ions for RNA input
  UI->seq_type->setCurrentText("RNA");
  p.clear();
  p = dialog_.getParam_();
  QVERIFY2(string(p.getValue("add_a-B_ions")) == "true", "Parameter 'add_a-B_ions' wasn't set correctly.");
  QVERIFY2(double(p.getValue("a-B_intensity")) == 1.23, "Parameter 'a-B_intensity' wasn't set correctly.");

  QVERIFY2(string(p.getValue("add_d_ions")) == "true", "Parameter 'add_d_ions' wasn't set correctly.");
  QVERIFY2(double(p.getValue("d_intensity")) == 1.23, "Parameter 'd_intensity' wasn't set correctly.");

  QVERIFY2(string(p.getValue("add_w_ions")) == "true", "Parameter 'add_w_ions' wasn't set correctly.");
  QVERIFY2(double(p.getValue("w_intensity")) == 1.23, "Parameter 'w_intensity' wasn't set correctly.");

  // Metabolite input doesn't have any exclusive parameters that need to be tested.
}

void TestTSGDialog::testSpectrumCalculation()
{
  // Spectrum generation settings don't really matter, since the dialog test shouldn't test the generators.
  // Hence, it is only checked if a spectrum was produced.
  
  UI->seq_type->setCurrentText("Peptide");
  QTest::qWait(DELAY);
  UI->seq_input->setText("PEPTIDE"); // set peptide sequence to 'PEPTIDE'
  QTest::qWait(DELAY);
  QTest::mouseClick(UI->dialog_buttons->button(QDialogButtonBox::Ok), Qt::LeftButton);
  QTest::qWait(DELAY);
  MSSpectrum pep_spec = dialog_.getSpectrum();
  QVERIFY2(!pep_spec.empty(), "Peptide input didn't produce a spectrum.");
  
  UI->seq_type->setCurrentText("RNA");
  QTest::qWait(DELAY);
  UI->seq_input->setText("AGUCCG");
  QTest::qWait(DELAY);
  QTest::mouseClick(UI->dialog_buttons->button(QDialogButtonBox::Ok), Qt::LeftButton);
  QTest::qWait(DELAY);
  MSSpectrum rna_spec = dialog_.getSpectrum();
  QVERIFY2(!rna_spec.empty(), "RNA input didn't produce a spectrum.");
  
  UI->seq_type->setCurrentText("Metabolite");
  QTest::qWait(DELAY);
  UI->seq_input->setText("C100H70N2O6");
  QTest::qWait(DELAY);
  QTest::mouseClick(UI->dialog_buttons->button(QDialogButtonBox::Ok), Qt::LeftButton);
  QTest::qWait(DELAY);
  MSSpectrum meta_spec = dialog_.getSpectrum();
  QVERIFY2(!meta_spec.empty(), "Metabolite input didn't produce a spectrum.");
}

void TestTSGDialog::testErrors()
{
  UI->seq_type->setCurrentText("Peptide");
  testMessageBoxes_();

  UI->seq_type->setCurrentText("RNA");
  testMessageBoxes_();

  UI->seq_type->setCurrentText("Metabolite");
  testMessageBoxes_();
}

// expands to a simple main() method that runs all the private slots
QTEST_MAIN(TestTSGDialog)

// Sadly this doesn't work and yields a seq fault, when trying to access any member of the dropDownList.
// That's why the combo box ('seq_type') is set manually and therefore not actually tested in this test.
// I was unable to fix this, but maybe someone else will get this to work..
// sources:
// https://gist.github.com/peteristhegreat/cbd8eaa0e565d0b82dbfb5c7fdc61c8d
// https://vicrucann.github.io/tutorials/qttest-signals-qtreewidget/
// https://stackoverflow.com/questions/69283103/qts-qtest-doesnt-select-an-item-in-a-drop-down-list-with-a-click
//
// #include <qcombobox.h>
// #include <qlistview.h>
//
// void TestTSGDialog::clickDropDown_(int row, QComboBox* comboBox)
//{
//  QListView* dropDownList = comboBox->findChild<QListView*>();
//  QModelIndex foundIndex {dropDownList->model()->index(row, 0)};
//
//  QRect foundDropDownItem = dropDownList->visualRect(foundIndex);
//  QPoint foundDropDownItemPosition = foundDropDownItem.center();
//
//  QWidget* activeWidget = dropDownList->viewport();
//  QTest::mouseClick(activeWidget, Qt::LeftButton, Qt::NoModifier, foundDropDownItemPosition);
//  QTest::qWait(DELAY);
//}

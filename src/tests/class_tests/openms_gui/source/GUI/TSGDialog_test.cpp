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
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <TSGDialog_test.h>
#include <ui_TheoreticalSpectrumGenerationDialog.h>

#include <qcombobox.h>
#include <qlistview.h>

using namespace OpenMS;
using namespace std;

constexpr int DELAY {1000};

// does the checkbox exist for <'Peptide','RNA'>
const map<Checkbox, pair<bool, bool>> intensity_ion_exists {
  {Checkbox::A_Ions, {1, 1}},
  {Checkbox::A_b_Ions, {0, 1}},
  {Checkbox::B_Ions, {1, 1}},
  {Checkbox::C_Ions, {1, 1}},
  {Checkbox::D_Ions, {0, 1}},
  {Checkbox::W_Ions, {0, 1}},
  {Checkbox::X_Ions, {1, 1}},
  {Checkbox::Y_Ions, {1, 1}},
  {Checkbox::Z_Ions, {1, 1}},
  {Checkbox::Precursor, {1, 1}},
  {Checkbox::Neutral_losses, {1, 0}},
  {Checkbox::Abundant_Immonium_Ions, {1, 0}}
};

// Sadly this doesn't work and yields a seq fault, when trying to access the dropDownList.
// That's why the spinbox is set manually and therefore not actually tested here.
// I was unable to fix this, but maybe someone else will get this to work..
// sources:
// https://gist.github.com/peteristhegreat/cbd8eaa0e565d0b82dbfb5c7fdc61c8d
// https://vicrucann.github.io/tutorials/qttest-signals-qtreewidget/
void clickDropDown(int row, QComboBox* comboBox)
{
  QListView* dropDownList = comboBox->findChild<QListView*>();
  QModelIndex foundIndex {dropDownList->model()->index(row, 0)};

  QRect foundDropDownItem = dropDownList->visualRect(foundIndex);
  QPoint foundDropDownItemPosition = foundDropDownItem.center();

  QWidget* activeWidget = dropDownList->viewport();
  QTest::mouseClick(activeWidget, Qt::LeftButton, Qt::NoModifier, foundDropDownItemPosition);
  QTest::qWait(DELAY); // waits 1 second
}

void TestTSGDialog::clickIsotopeModel_()
{
  // QTest::mouseClick needs the exact position of the interactable part of the button
  UI->model_none->click();
  QTest::qWait(DELAY);
  QVERIFY(!(UI->max_iso_spinbox->isEnabled()));
  QVERIFY(!(UI->max_iso_label->isEnabled()));
  QVERIFY(!(UI->max_iso_prob_spinbox->isEnabled()));
  QVERIFY(!(UI->max_iso_prob_label->isEnabled()));

  UI->model_coarse->click();
  QTest::qWait(DELAY);
  QVERIFY(UI->max_iso_spinbox->isEnabled());
  QVERIFY(UI->max_iso_label->isEnabled());
  QVERIFY(!(UI->max_iso_prob_spinbox->isEnabled()));
  QVERIFY(!(UI->max_iso_prob_label->isEnabled()));

  UI->model_fine->click();
  QTest::qWait(DELAY);
  QVERIFY(!(UI->max_iso_spinbox->isEnabled()));
  QVERIFY(!(UI->max_iso_label->isEnabled()));
  QVERIFY(UI->max_iso_prob_spinbox->isEnabled());
  QVERIFY(UI->max_iso_prob_label->isEnabled());
}

void TestTSGDialog::testConstruction()
{
  // editable/interactable GUI parts
  QVERIFY2(UI->seq_type, "Sequence selection combo box not created.");
  QVERIFY2(UI->seq_input, "Sequence input line edit not created.");
  QVERIFY2(UI->charge_spinbox, "Charge spin box not created.");
  QVERIFY2(UI->max_iso_spinbox, "Max. isotope model spin box not created.");
  QVERIFY2(UI->max_iso_prob_spinbox, "Max. isotope probability spin box not created.");
  QVERIFY2(UI->list_widget, "Ion list widget not created.");
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
  QVERIFY2(UI->button_box, "Buttonbox not created.");

  // labels
  QVERIFY2(UI->enter_seq_label, "'Enter sequence' label not created.");
  QVERIFY2(UI->charge_label, "'Charge' label not created.");
  QVERIFY2(UI->generate_label, "'Generate' label not created.");
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

void TestTSGDialog::testSpectrumCalculation()
{
}

void TestTSGDialog::testGui()
{
  const map<Checkbox, QDoubleSpinBox*> checkbox_to_intensity_ {{Checkbox::A_Ions, UI->a_intensity},
                                                               {Checkbox::A_b_Ions, UI->a_b_intensity},
                                                               {Checkbox::B_Ions, UI->b_intensity},
                                                               {Checkbox::C_Ions, UI->c_intensity},
                                                               {Checkbox::D_Ions, UI->d_intensity},
                                                               {Checkbox::W_Ions, UI->w_intensity},
                                                               {Checkbox::X_Ions, UI->x_intensity},
                                                               {Checkbox::Y_Ions, UI->y_intensity},
                                                               {Checkbox::Z_Ions, UI->z_intensity},
                                                               {Checkbox::Precursor, nullptr},
                                                               {Checkbox::Neutral_losses, nullptr}, // UI->rel_loss_intensity is a normal spin box
                                                               {Checkbox::Abundant_Immonium_Ions, nullptr}};

  dialog_.show();

  //////////////////////////////////////////////////////
  //                     PEPTIDE                      //
  //////////////////////////////////////////////////////
  UI->seq_type->setCurrentText("Peptide");
  QTest::qWait(DELAY);

  // isotope model
  QVERIFY2(!UI->isotope_model->isHidden(), "Isotope model hidden for 'Peptide' setting.");
  clickIsotopeModel_();

  // ion types and intensities
  for (const Checkbox& c : check_box_names)
  {
    // get the item
    QListWidgetItem* item = UI->list_widget->item(int(c));
    QVERIFY(item);

    // get intensity spin box corresponding to current check box
    QDoubleSpinBox* spin = checkbox_to_intensity_.at(c);
    
    if(intensity_ion_exists.at(c).first)
    {
      // check state before clicking
      Qt::CheckState prev = item->checkState();

      // get the rectangular coordinates of the item
      QRect rect = UI->list_widget->visualItemRect(item);

      // imitate the click on checkbox c
      QTest::mouseClick(UI->list_widget->viewport(), Qt::LeftButton, 0, rect.center());
      QTest::qWait(DELAY);

      // verfiy the check state changed
      QVERIFY(prev != item->checkState());

      if (spin == nullptr) continue;

      // simulate keyboard input
      spin->clear();
      QTest::keyClicks(spin, "2");
      QTest::qWait(DELAY);

      QVERIFY(2.0 == spin->value());
    }
    else
    {
      // if ion type isn't supported, check if the ion and its intensity are hidden
      QVERIFY(item->isHidden());
      if (spin == nullptr) continue;
      QVERIFY(spin->isHidden());
    }
  }
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

  // isotope model
  QVERIFY2(UI->isotope_model->isHidden(), "Isotope model not hidden for 'Peptide' setting.");
  QVERIFY(UI->max_iso_spinbox->isHidden());
  QVERIFY(UI->max_iso_label->isHidden());
  QVERIFY(UI->max_iso_prob_spinbox->isHidden());
  QVERIFY(UI->max_iso_prob_label->isHidden());

  // ion types and intensities
  for (const Checkbox& c : check_box_names)
  {
    // get the item
    QListWidgetItem* item = UI->list_widget->item(int(c));
    QVERIFY(item);

    // get intensity spin box corresponding to current check box
    QDoubleSpinBox* spin = checkbox_to_intensity_.at(c);

    if (intensity_ion_exists.at(c).second)
    {
      // check state before clicking
      Qt::CheckState prev = item->checkState();

      // get the rectangular coordinates of the item
      QRect rect = UI->list_widget->visualItemRect(item);

      // imitate the click on checkbox c
      QTest::mouseClick(UI->list_widget->viewport(), Qt::LeftButton, 0, rect.center());
      QTest::qWait(DELAY);

      // verfiy the check state changed
      QVERIFY(prev != item->checkState());

      if (spin == nullptr)
        continue;

      // simulate keyboard input
      spin->clear();
      QTest::keyClicks(spin, "2");
      QTest::qWait(DELAY);

      QVERIFY(2.0 == spin->value());
    }
    else
    {
      // if ion type isn't supported, check if the ion and its intensity are hidden
      QVERIFY(item->isHidden());
      if (spin == nullptr) continue;
      QVERIFY(spin->isHidden());
    }
  }
  // check relative loss intensity manually
  QVERIFY(UI->rel_loss_intensity->isHidden());
  QVERIFY(UI->rel_loss_label->isHidden());

}

// expands to a simple main() method that runs all the private slots (test functions)
QTEST_MAIN(TestTSGDialog)
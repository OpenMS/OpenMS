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

#define UI dialog_.ui_

using namespace OpenMS;

constexpr int DELAY {1000};

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
  QTest::mouseClick(UI->model_none, Qt::LeftButton);
  QVERIFY(!(UI->max_iso_spinbox->isEnabled()));
  QVERIFY(!(UI->max_iso_label->isEnabled()));
  QVERIFY(!(UI->max_iso_prob_spinbox->isEnabled()));
  QVERIFY(!(UI->max_iso_prob_label->isEnabled()));

  QTest::mouseClick(UI->model_coarse, Qt::LeftButton);
  QVERIFY(UI->max_iso_spinbox->isEnabled());
  QVERIFY(UI->max_iso_label->isEnabled());
  QVERIFY(!(UI->max_iso_prob_spinbox->isEnabled()));
  QVERIFY(!(UI->max_iso_prob_label->isEnabled()));

  QTest::mouseClick(UI->model_fine, Qt::LeftButton);
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
  UI->seq_type->setCurrentText("Peptide");
  clickIsotopeModel_();

  UI->seq_type->setCurrentText("RNA");

  UI->seq_type->setCurrentText("Metabolite");
  clickIsotopeModel_();
}

// expands to a simple main() method that runs all the test functions
QTEST_MAIN(TestTSGDialog)
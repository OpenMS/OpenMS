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

constexpr int DELAY {1000};

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

void TestTSGDialog::testConstruction()
{
  // editable/interactable GUI parts
  QVERIFY(dialog_.ui_->seq_type, "Sequence selection combo box not created.");
  QVERIFY(dialog_.ui_->seq_input, "Sequence input line edit not created.");
  QVERIFY(dialog_.ui_->charge_spinbox, "Charge spin box not created.");
  QVERIFY(dialog_.ui_->max_iso_spinbox, "Max. isotope model spin box not created.");
  QVERIFY(dialog_.ui_->max_iso_prob_spinbox, "Max. isotope probability spin box not created.");
  QVERIFY(dialog_.ui_->list_widget, "Ion list widget not created.");
  QVERIFY(dialog_.ui_->a_intensity, "A ion intensity spin box not created.");
  QVERIFY(dialog_.ui_->a_b_intensity, "A-b ion intensity spin box not created.");
  QVERIFY(dialog_.ui_->b_intensity, "B ion intensity spin box not created.");
  QVERIFY(dialog_.ui_->c_intensity, "C ion intensity spin box not created.");
  QVERIFY(dialog_.ui_->d_intensity, "D ion intensity spin box not created.");
  QVERIFY(dialog_.ui_->w_intensity, "W ion intensity spin box not created.");
  QVERIFY(dialog_.ui_->x_intensity, "X ion intensity spin box not created.");
  QVERIFY(dialog_.ui_->y_intensity, "Y ion intensity spin box not created.");
  QVERIFY(dialog_.ui_->z_intensity, "Z ion intensity spin box not created.");
  QVERIFY(dialog_.ui_->rel_loss_intensity, "Relative loss intensity spin box not created.");
  QVERIFY(dialog_.ui_->button_box, "Buttonbox not created.");

  // labels
  QVERIFY(dialog_.ui_->enter_seq_label, "'Enter sequence' label not created.");
  QVERIFY(dialog_.ui_->charge_label, "'Charge' label not created.");
  QVERIFY(dialog_.ui_->generate_label, "'Generate' label not created.");
  QVERIFY(dialog_.ui_->max_iso_label, "'Max. Isotope' label not created.");
  QVERIFY(dialog_.ui_->max_iso_prob_label, "'Max. Isotope Probability in %' label not created.");
  QVERIFY(dialog_.ui_->a_label, "'A-ions' label not created.");
  QVERIFY(dialog_.ui_->a_b_label, "'A-b-ions' label not created.");
  QVERIFY(dialog_.ui_->b_label, "'B-ions' label not created.");
  QVERIFY(dialog_.ui_->c_label, "'C-ions' label not created.");
  QVERIFY(dialog_.ui_->d_label, "'D-ions' label not created.");
  QVERIFY(dialog_.ui_->w_label, "'W-ions' label not created.");
  QVERIFY(dialog_.ui_->x_label, "'X-ions' label not created.");
  QVERIFY(dialog_.ui_->y_label, "'Y-ions' label not created.");
  QVERIFY(dialog_.ui_->z_label, "'Z-ions' label not created.");
  QVERIFY(dialog_.ui_->rel_loss_label, "'Relative loss in %' label not created.");
  
  // group boxes
  QVERIFY(dialog_.ui_->isotope_model, "Isotope model group box not created.");
  QVERIFY(dialog_.ui_->intensities, "Intensity group box not created.");
}

void TestTSGDialog::testGui()
{
  //QVERIFY(dialog_.ui_->button_box, "Buttonbox not created.");
  clickDropDown(0, dialog_.ui_->seq_type); // select 'Peptide'
  clickDropDown(1, dialog_.ui_->seq_type); // select 'RNA'
}

// expands to a simple main() method that runs all the test functions
QTEST_MAIN(TestTSGDialog)
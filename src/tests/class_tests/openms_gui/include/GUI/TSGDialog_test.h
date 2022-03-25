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

#include <QtTest/QtTest>
#include <QtGui>
#include <qspinbox.h>

#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>
#include <ui_TheoreticalSpectrumGenerationDialog.h>

#define UI dialog_.ui_

namespace OpenMS
{
  class TestTSGDialog : public QObject
  {
    Q_OBJECT

    public:
    TestTSGDialog() : dialog_() {}

    ~TestTSGDialog()
    {
      dialog_.destroy();
    }

    private slots:
      void testConstruction();
      
      void testGui();

      void testParameterImport();

      void testSpectrumCalculation();

    private:
      template<typename T> // template for QSpinBox and QDoubleSpinBox
      void testSpinBox_(T* box, std::string str_value = "2");

      void testIonsIntensities_(bool peptide_input);

      void testSequenceInput_(QString input);

      void testIsotopeModel_();

      TheoreticalSpectrumGenerationDialog dialog_;

      // For each check box get its intensity spin box.
      // Order is important here!
      // To access the right entry for each check box
      // use int(TheoreticalSpectrumGenerationDialog::CheckBox).
      const std::vector<QDoubleSpinBox*> check_box_to_intensity_ {
        UI->a_intensity,        // A-Ion
        UI->a_b_intensity,      // a-B-Ion
        UI->b_intensity,        // B-Ion
        UI->c_intensity,        // C-Ion
        UI->d_intensity,        // D-Ion
        UI->w_intensity,        // W-Ion
        UI->x_intensity,        // X-Ion
        UI->y_intensity,        // Y-Ion
        UI->z_intensity,        // Z-Ion
        nullptr,                // Precursor
        nullptr,                // Neutral losses: UI->rel_loss_intensity is a normal spin box
        nullptr                 // Abundant Immonium Ions
      };

      // To check if the check box is enabled for <'Peptide','RNA'>
      // Order is important here!
      // To access the right entry for each check box
      // use int(TheoreticalSpectrumGenerationDialog::CheckBox).
      const std::vector<std::pair<bool, bool>> intensity_ion_exists {
        {1, 1},        // A-Ion
        {0, 1},        // a-B-Ion
        {1, 1},        // B-Ion
        {1, 1},        // C-Ion
        {0, 1},        // D-Ion
        {0, 1},        // W-Ion
        {1, 1},        // X-Ion
        {1, 1},        // Y-Ion
        {1, 1},        // Z-Ion
        {1, 1},        // Precursor
        {1, 0},        // Neutral losses
        {1, 0}         // Abundant Immonium Ions
       };
  };
}
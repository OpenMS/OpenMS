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

      // for each check box map its intensity spin box
      const std::map<Checkbox, QDoubleSpinBox*> checkbox_to_intensity_ {
        {Checkbox::A_Ions, UI->a_intensity},
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
        {Checkbox::Abundant_Immonium_Ions, nullptr}
      };

      // is the check box enabled for <'Peptide','RNA'>
      const std::map<Checkbox, std::pair<bool, bool>> intensity_ion_exists {
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
  };
}
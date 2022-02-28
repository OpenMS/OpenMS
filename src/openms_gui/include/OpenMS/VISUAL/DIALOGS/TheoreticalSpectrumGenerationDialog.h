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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#pragma once

#include <QtWidgets/QDialog>
class QListWidgetItem;
class QListWidget;
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace Ui
{
  class TheoreticalSpectrumGenerationDialogTemplate;
}

namespace OpenMS
{
  /**
      @brief Dialog which allows to enter an AA sequence and generates a theoretical spectrum for it.

      @ingroup Dialogs
  */
  class TheoreticalSpectrumGenerationDialog :
    public QDialog
  {
    Q_OBJECT

public:

    // Note: If an additional checkbox is added all of the following three objects have to be edited!
    // 
    // 
    // enum to get ion checkbox index (ordering is important here!)
    enum class Checkbox
    {
      A_Ions,
      B_Ions,
      C_Ions,
      X_Ions,
      Y_Ions,
      Z_Ions,
      Precursor,
      Neutral_losses,
      Isotope_cluster,
      Abundant_Immonium_Ions
    };

    // vector of all enum entries for iteration
    const std::vector<Checkbox> check_box_names {Checkbox::A_Ions,
                                                 Checkbox::B_Ions,
                                                 Checkbox::C_Ions,
                                                 Checkbox::X_Ions,
                                                 Checkbox::Y_Ions,
                                                 Checkbox::Z_Ions,
                                                 Checkbox::Precursor,
                                                 Checkbox::Neutral_losses,
                                                 Checkbox::Isotope_cluster,
                                                 Checkbox::Abundant_Immonium_Ions};

    // map from checkbox (index) to corresponding parameter with description
    const std::map<Checkbox, std::pair<String, String>> checkbox_to_param {
      {Checkbox::A_Ions, {"add_a_ions", "Add peaks of a-ions to the spectrum"}},
      {Checkbox::B_Ions, {"add_b_ions", "Add peaks of b-ions to the spectrum"}},
      {Checkbox::C_Ions, {"add_c_ions", "Add peaks of c-ions to the spectrum"}},
      {Checkbox::X_Ions, {"add_x_ions", "Add peaks of x-ions to the spectrum"}},
      {Checkbox::Y_Ions, {"add_y_ions", "Add peaks of y-ions to the spectrum"}},
      {Checkbox::Z_Ions, {"add_z_ions", "Add peaks of z-ions to the spectrum"}},
      {Checkbox::Precursor, {"add_precursor_peaks", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes"}},
      {Checkbox::Neutral_losses, {"add_losses", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered"}},
      {Checkbox::Isotope_cluster, {"add_isotopes", "If set to 1 isotope peaks of the product ion peaks are added"}},
      {Checkbox::Abundant_Immonium_Ions, {"add_abundant_immonium_ions", "Add most abundant immonium ions"}}
    };

    /// Constructor
    TheoreticalSpectrumGenerationDialog();
    /// Destructor
    ~TheoreticalSpectrumGenerationDialog() override;

    String getSequence() const;

    Param getParam() const;

protected slots:

    void itemChanged(QListWidgetItem * item);

protected:

private:
    Ui::TheoreticalSpectrumGenerationDialogTemplate* ui_;
  };

}

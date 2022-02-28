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

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>
#include <ui_TheoreticalSpectrumGenerationDialog.h>

#include <OpenMS/DATASTRUCTURES/String.h>


namespace OpenMS
{
  TheoreticalSpectrumGenerationDialog::TheoreticalSpectrumGenerationDialog()
    : ui_(new Ui::TheoreticalSpectrumGenerationDialogTemplate)
  {
    ui_->setupUi(this);

    connect(ui_->list_widget, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(itemChanged(QListWidgetItem *)));

    // select b- and y-ions as residue types by default
    for (Checkbox c : check_box_names)
    {
      if (c == Checkbox::B_Ions || c == Checkbox::Y_Ions)
      {
        ui_->list_widget->item(int(c))->setCheckState(Qt::Checked);
        continue;
      }
      ui_->list_widget->item(int(c))->setCheckState(Qt::Unchecked);
    }
  }

  TheoreticalSpectrumGenerationDialog::~TheoreticalSpectrumGenerationDialog()
  {
    delete ui_;
  }

  String TheoreticalSpectrumGenerationDialog::getSequence() const
  {
    return ui_->line_edit->text();
  }

  Param TheoreticalSpectrumGenerationDialog::getParam() const
  {
    Param p;
    p.setValue("charge", ui_->spin_box->value());

    // add checkboxes to parameters, i.e. ion types
    for (Checkbox c : check_box_names)
    {
      bool status = (ui_->list_widget->item(int(c))->checkState() == Qt::Checked);
      String status_str = status ? "true" : "false";
      p.setValue(checkbox_to_param.at(c).first, status_str, checkbox_to_param.at(c).second);
    }

    // add intensities
    Size max_iso_count = (Size)ui_->max_iso_spinbox->value();
    p.setValue("max_isotope", max_iso_count, "Number of isotopic peaks");
    p.setValue("a_intensity", ui_->a_intensity->value(), "Intensity of the a-ions");
    p.setValue("b_intensity", ui_->b_intensity->value(), "Intensity of the b-ions");
    p.setValue("c_intensity", ui_->c_intensity->value(), "Intensity of the c-ions");
    p.setValue("x_intensity", ui_->x_intensity->value(), "Intensity of the x-ions");
    p.setValue("y_intensity", ui_->y_intensity->value(), "Intensity of the y-ions");
    p.setValue("z_intensity", ui_->z_intensity->value(), "Intensity of the z-ions");
    double rel_loss_int = (double)(ui_->rel_loss_intensity->value()) / 100.0;
    p.setValue("relative_loss_intensity", rel_loss_int, "Intensity of loss ions, in relation to the intact ion intensity");

    return p;
  }

  void TheoreticalSpectrumGenerationDialog::itemChanged(QListWidgetItem * item)
  {
    if (item->text() == "Isotope clusters")
    {
      if (item->checkState() == Qt::Checked)
      {
        ui_->max_iso_label->setEnabled(true);
        ui_->max_iso_spinbox->setEnabled(true);
      }
      else
      {
        ui_->max_iso_label->setEnabled(false);
        ui_->max_iso_spinbox->setEnabled(false);
      }
    }
  }

} // namespace

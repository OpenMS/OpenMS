// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
    ui_->list_widget->item(0)->setCheckState(Qt::Unchecked);
    ui_->list_widget->item(1)->setCheckState(Qt::Checked);
    ui_->list_widget->item(2)->setCheckState(Qt::Unchecked);
    ui_->list_widget->item(3)->setCheckState(Qt::Unchecked);
    ui_->list_widget->item(4)->setCheckState(Qt::Checked);
    ui_->list_widget->item(5)->setCheckState(Qt::Unchecked);
    ui_->list_widget->item(6)->setCheckState(Qt::Unchecked);
    ui_->list_widget->item(7)->setCheckState(Qt::Unchecked);
    ui_->list_widget->item(8)->setCheckState(Qt::Unchecked);
    ui_->list_widget->item(9)->setCheckState(Qt::Unchecked);
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
    bool losses = (ui_->list_widget->item(7)->checkState() == Qt::Checked); // "Neutral losses"
    String losses_str = losses ? "true" : "false";
    p.setValue("add_losses", losses_str, "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");

    bool isotopes = (ui_->list_widget->item(8)->checkState() == Qt::Checked); // "Isotope clusters"
    String iso_str = isotopes ? "true" : "false";
    p.setValue("add_isotopes", iso_str, "If set to 1 isotope peaks of the product ion peaks are added");

    bool abundant_immonium_ions = (ui_->list_widget->item(9)->checkState() == Qt::Checked); // "abundant immonium-ions"
    String abundant_immonium_ions_str = abundant_immonium_ions ? "true" : "false";
    p.setValue("add_abundant_immonium_ions", abundant_immonium_ions_str, "Add most abundant immonium ions");

    bool precursor_ions = (ui_->list_widget->item(6)->checkState() == Qt::Checked); // "add precursor ions"
    String precursor_ions_str = precursor_ions ? "true" : "false";
    p.setValue("add_precursor_peaks", precursor_ions_str, "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");

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
    
    p.setValue("has_A", "false");
    p.setValue("has_B", "false");
    p.setValue("has_C", "false");
    p.setValue("has_X", "false");
    p.setValue("has_Y", "false");
    p.setValue("has_Z", "false");
    p.setValue("has_Precursor", "false");
    p.setValue("has_abundantImmoniumIons", "false");

    if (ui_->list_widget->item(0)->checkState() == Qt::Checked) // "A-ions"
    {
      p.setValue("has_A", "true");
    }
    if (ui_->list_widget->item(1)->checkState() == Qt::Checked) // "B-ions"
    {
      p.setValue("has_B", "true");
    }
    if (ui_->list_widget->item(2)->checkState() == Qt::Checked) // "C-ions"
    {
      p.setValue("has_C", "true");
    }
    if (ui_->list_widget->item(3)->checkState() == Qt::Checked) // "X-ions"
    {
      p.setValue("has_X", "true");
    }
    if (ui_->list_widget->item(4)->checkState() == Qt::Checked) // "Y-ions"
    {
      p.setValue("has_Y", "true");
    }
    if (ui_->list_widget->item(5)->checkState() == Qt::Checked) // "Z-ions"
    {
      p.setValue("has_Z", "true");
    }
    if (ui_->list_widget->item(6)->checkState() == Qt::Checked) // "Precursor"
    {
      p.setValue("has_Precursor", "true");
    }
    if (ui_->list_widget->item(9)->checkState() == Qt::Checked) // "abundant Immonium-ions"
    {
      p.setValue("has_abundantImmoniumIons", "true");
    }
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

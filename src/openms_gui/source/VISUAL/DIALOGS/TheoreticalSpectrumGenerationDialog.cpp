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

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/String.h>

// Qt includes
#include <QtWidgets/QMessageBox>


namespace OpenMS
{
  TheoreticalSpectrumGenerationDialog::TheoreticalSpectrumGenerationDialog()
    : ui_(new Ui::TheoreticalSpectrumGenerationDialogTemplate)
  {
    ui_->setupUi(this);

    // if dialog is accepted, try generating a spectrum, only close dialog on success
    connect(ui_->button_box, SIGNAL(accepted()), this, SLOT(calculateSpectrum()));

    // signals for changing isotope model interface
    // there has to be a smarter way of doing this
    connect(ui_->model_none, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));
    connect(ui_->model_coarse, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));
    connect(ui_->model_fine, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));

    // don't add any isotopes by default
    ui_->model_none->setChecked(true);

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

    // add checkboxes to parameters, i.e. ion types
    for (Checkbox c : check_box_names)
    {
      bool status = (ui_->list_widget->item(int(c))->checkState() == Qt::Checked);
      String status_str = status ? "true" : "false";
      p.setValue(checkbox_to_param.at(c).first, status_str, checkbox_to_param.at(c).second);
    }

    // charge
    p.setValue("charge", ui_->spin_box->value());

    // isotopes
    if (!ui_->model_none->isChecked()) // add isotopes if any other model than 'None' is chosen
    {
      bool coarse_model = ui_->model_coarse->isChecked();
      String model = coarse_model ? "coarse" : "fine";
      p.setValue("isotope_model", model , "Model to use for isotopic peaks ('none' means no isotopic peaks are added, 'coarse' adds isotopic peaks in unit mass distance, 'fine' uses the hyperfine isotopic generator to add accurate isotopic peaks. Note that adding isotopic peaks is very slow.");

      if (coarse_model)
      {
        Size max_iso_count = (Size)ui_->max_iso_spinbox->value();
        p.setValue("max_isotope", max_iso_count, "Defines the maximal isotopic peak which is added if 'isotope_model' is 'coarse'");
      }
      else
      {
        double max_iso_prob = (double)ui_->max_iso_prob_spinbox->value() / 100.;
        p.setValue("max_isotope_probability", max_iso_prob, "Defines the maximal isotopic probability to cover if 'isotope_model' is 'fine'");
      }
    }
    else // don't add isotopes
    {
      p.setValue("isotope_model", "none", "Model to use for isotopic peaks ('none' means no isotopic peaks are added, 'coarse' adds isotopic peaks in unit mass distance, 'fine' uses the hyperfine isotopic generator to add accurate isotopic peaks. Note that adding isotopic peaks is very slow.");
    }

    // add intensities
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

  MSSpectrum TheoreticalSpectrumGenerationDialog::getSpectrum() const
  {
    return spec_;
  }

  void TheoreticalSpectrumGenerationDialog::calculateSpectrum()
  {
    String seq_string(this->getSequence());
    if (seq_string == "")
    {
      QMessageBox::warning(this, "Error", "You must enter a peptide sequence!");
      return;
    }
    AASequence aa_sequence;
    try
    {
      aa_sequence = AASequence::fromString(seq_string);
    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::warning(this, "Error", QString("Spectrum generation failed! (") + e.what() + ")");
      return;
    }

    Param p = this->getParam();
    Int charge = p.getValue("charge");
    p.remove("charge"); // "charge" isn't a parameter of TheoreticalSpectrumGenerator

    p.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");

    TheoreticalSpectrumGenerator generator;
    generator.setParameters(p);

    try
    {
      generator.getSpectrum(spec_, aa_sequence, charge, charge);
    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::warning(this, "Error", QString("Spectrum generation failed! (") + e.what() + "). Please report this to the developers (specify what input you used)!");
      return;
    }

    if (spec_.empty())
    {
      QMessageBox::warning(this, "Error", QString("The generated spectrum was empty and will not be drawn!"));
      return;
    }

    this->accept();
  }

  void TheoreticalSpectrumGenerationDialog::modelChanged()
  {    
    if (ui_->model_none->isChecked())
    {
      ui_->max_iso_label->setEnabled(false);
      ui_->max_iso_spinbox->setEnabled(false);
      ui_->max_iso_prob_label->setEnabled(false);
      ui_->max_iso_prob_spinbox->setEnabled(false);
    }
    else if (ui_->model_coarse->isChecked())
    {
      ui_->max_iso_label->setEnabled(true);
      ui_->max_iso_spinbox->setEnabled(true);
      ui_->max_iso_prob_label->setEnabled(false);
      ui_->max_iso_prob_spinbox->setEnabled(false);
    }
    else if (ui_->model_fine->isChecked())
    {
      ui_->max_iso_label->setEnabled(false);
      ui_->max_iso_spinbox->setEnabled(false);
      ui_->max_iso_prob_label->setEnabled(true);
      ui_->max_iso_prob_spinbox->setEnabled(true);
    }

  }

} // namespace

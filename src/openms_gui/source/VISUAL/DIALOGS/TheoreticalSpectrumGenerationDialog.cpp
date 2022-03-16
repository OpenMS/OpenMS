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
// $Authors: Marc Sturm, Tom Waschischeck $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/String.h>

// Qt includes
#include <QtWidgets/QMessageBox>
#include <qflags.h>


namespace OpenMS
{
  // Note: If an additional checkbox is added all of the following three objects have to be edited!
  //
  //
  // enum to get ion checkbox index (Ordering has to be the same as in the ui!)
  enum class Checkbox
  {
    A_Ions,
    A_b_Ions,
    B_Ions,
    C_Ions,
    D_Ions,
    W_Ions,
    X_Ions,
    Y_Ions,
    Z_Ions,
    Precursor,
    Neutral_losses,
    Abundant_Immonium_Ions
  };

  // vector of all enum entries for iteration
  const std::vector<Checkbox> check_box_names {Checkbox::A_Ions, Checkbox::A_b_Ions, Checkbox::B_Ions, Checkbox::C_Ions,    Checkbox::D_Ions,         Checkbox::W_Ions,
                                               Checkbox::X_Ions, Checkbox::Y_Ions,   Checkbox::Z_Ions, Checkbox::Precursor, Checkbox::Neutral_losses, Checkbox::Abundant_Immonium_Ions};

  // map from checkbox (index) to corresponding parameter with description
  const std::map<Checkbox, std::pair<String, String>> checkbox_to_param {
    {Checkbox::A_Ions, {"add_a_ions", "Add peaks of a-ions to the spectrum"}},
    {Checkbox::A_b_Ions, {"add_a-B_ions", "Add peaks of a-B-ions to the spectrum (nucleotide sequences only)"}},
    {Checkbox::B_Ions, {"add_b_ions", "Add peaks of b-ions to the spectrum"}},
    {Checkbox::C_Ions, {"add_c_ions", "Add peaks of c-ions to the spectrum"}},
    {Checkbox::D_Ions, {"add_d_ions", "Add peaks of d-ions to the spectrum (nucleotide sequences only)"}},
    {Checkbox::W_Ions, {"add_w_ions", "Add peaks of w-ions to the spectrum (nucleotide sequences only)"}},
    {Checkbox::X_Ions, {"add_x_ions", "Add peaks of x-ions to the spectrum"}},
    {Checkbox::Y_Ions, {"add_y_ions", "Add peaks of y-ions to the spectrum"}},
    {Checkbox::Z_Ions, {"add_z_ions", "Add peaks of z-ions to the spectrum"}},
    {Checkbox::Precursor, {"add_precursor_peaks", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes"}},
    {Checkbox::Neutral_losses, {"add_losses", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered (peptide sequences only)"}},
    {Checkbox::Abundant_Immonium_Ions, {"add_abundant_immonium_ions", "Add most abundant immonium ions (peptide sequences only)"}}};

  // specific checkboxes (TheoreticalSpectrumGenerator (peptide) vs. NucleicAcidSpectrumGenerator (rna))
  const std::vector<Checkbox> rna_specific_ions {Checkbox::A_b_Ions, Checkbox::D_Ions, Checkbox::W_Ions};
  const std::vector<Checkbox> peptide_specific_ions {Checkbox::Neutral_losses, Checkbox::Abundant_Immonium_Ions};

  TheoreticalSpectrumGenerationDialog::TheoreticalSpectrumGenerationDialog()
  {
    ui_.setupUi(this);

    // if dialog is accepted, try generating a spectrum, only close dialog on success
    connect(ui_.button_box, &QDialogButtonBox::accepted, this, &TheoreticalSpectrumGenerationDialog::calculateSpectrum);

    // signals for changing isotope model interface
    connect(ui_.model_none, &QRadioButton::toggled, this, &TheoreticalSpectrumGenerationDialog::modelChanged);
    connect(ui_.model_coarse, &QRadioButton::toggled, this, &TheoreticalSpectrumGenerationDialog::modelChanged);
    connect(ui_.model_fine, &QRadioButton::toggled, this, &TheoreticalSpectrumGenerationDialog::modelChanged);

    // don't add any isotopes by default and update interface
    ui_.model_none->setChecked(true);
    modelChanged();

    // signal for changing interface depending on sequence type
    connect(ui_.seq_type, &QComboBox::currentTextChanged, this, &TheoreticalSpectrumGenerationDialog::seqTypeSwitch);

    // select peptide sequence by default and update interface
    ui_.seq_type->setCurrentText("Peptide");
    seqTypeSwitch();

    // select b- and y-ions as residue types by default
    for (const Checkbox& c : check_box_names)
    {
      if (c == Checkbox::B_Ions || c == Checkbox::Y_Ions)
      {
        ui_.list_widget->item(int(c))->setCheckState(Qt::Checked);
        continue;
      }
      ui_.list_widget->item(int(c))->setCheckState(Qt::Unchecked);
    }

    // automatic layout
    // disables manual resizing from the user
    layout()->setSizeConstraint(QLayout::SetFixedSize);
  }

  TheoreticalSpectrumGenerationDialog::~TheoreticalSpectrumGenerationDialog() {}

  String TheoreticalSpectrumGenerationDialog::getSequence() const
  {
    return ui_.line_edit->text();
  }

  Param TheoreticalSpectrumGenerationDialog::getParam() const
  {
    Param p;

    bool peptide_input = ui_.seq_type->currentText() == "Peptide";

    // add checkboxes to parameters, i.e. ion types
    for (const Checkbox& c : check_box_names)
    {
      // for peptide input skip rna specific ions
      if (peptide_input && (std::find(rna_specific_ions.begin(), rna_specific_ions.end(), c) != rna_specific_ions.end())) continue;

      // for rna input skip peptide specific ions
      if (!peptide_input && (std::find(peptide_specific_ions.begin(), peptide_specific_ions.end(), c) != peptide_specific_ions.end())) continue;

      bool status = (ui_.list_widget->item(int(c))->checkState() == Qt::Checked);
      String status_str = status ? "true" : "false";
      p.setValue(checkbox_to_param.at(c).first, status_str, checkbox_to_param.at(c).second);
    }

    // charge
    p.setValue("charge", ui_.spin_box->value());

    // add intensities
    p.setValue("a_intensity", ui_.a_intensity->value(), "Intensity of the a-ions");
    p.setValue("b_intensity", ui_.b_intensity->value(), "Intensity of the b-ions");
    p.setValue("c_intensity", ui_.c_intensity->value(), "Intensity of the c-ions");
    p.setValue("x_intensity", ui_.x_intensity->value(), "Intensity of the x-ions");
    p.setValue("y_intensity", ui_.y_intensity->value(), "Intensity of the y-ions");
    p.setValue("z_intensity", ui_.z_intensity->value(), "Intensity of the z-ions");

    if (peptide_input) // peptide specific settings (TheoreticalSpectrumGenerator)
    {
      // isotopes
      if (!ui_.model_none->isChecked()) // add isotopes if any other model than 'None' is chosen
      {
        bool coarse_model = ui_.model_coarse->isChecked();
        String model = coarse_model ? "coarse" : "fine";
        p.setValue("isotope_model", model,
                   "Model to use for isotopic peaks ('none' means no isotopic peaks are added, 'coarse' adds isotopic peaks in unit mass distance, 'fine' uses the hyperfine isotopic generator to add "
                   "accurate isotopic peaks. Note that adding isotopic peaks is very slow.");

        if (coarse_model)
        {
          Size max_iso_count = (Size)ui_.max_iso_spinbox->value();
          p.setValue("max_isotope", max_iso_count, "Defines the maximal isotopic peak which is added if 'isotope_model' is 'coarse'");
        }
        else
        {
          double max_iso_prob = (double)ui_.max_iso_prob_spinbox->value() / 100.;
          p.setValue("max_isotope_probability", max_iso_prob, "Defines the maximal isotopic probability to cover if 'isotope_model' is 'fine'");
        }
      }
      else // don't add isotopes
      {
        p.setValue("isotope_model", "none",
                   "Model to use for isotopic peaks ('none' means no isotopic peaks are added, 'coarse' adds isotopic peaks in unit mass distance, 'fine' uses the hyperfine isotopic generator to add "
                   "accurate isotopic peaks. Note that adding isotopic peaks is very slow.");
      }

      // loss intensity
      double rel_loss_int = (double)(ui_.rel_loss_intensity->value()) / 100.0;
      p.setValue("relative_loss_intensity", rel_loss_int, "Intensity of loss ions, in relation to the intact ion intensity");
    }
    else // rna specific settings (NucleicAcidSpectrumGenerator)
    {
      // specific ion intensities
      p.setValue("a-B_intensity", ui_.a_b_intensity->value(), "Intensity of the a-B-ions");
      p.setValue("d_intensity", ui_.d_intensity->value(), "Intensity of the d-ions");
      p.setValue("w_intensity", ui_.w_intensity->value(), "Intensity of the w-ions");
    }

    return p;
  }

  const MSSpectrum& TheoreticalSpectrumGenerationDialog::getSpectrum() const
  {
    return spec_;
  }

  void TheoreticalSpectrumGenerationDialog::calculateSpectrum()
  {
    bool peptide_input = ui_.seq_type->currentText() == "Peptide";

    String seq_string(this->getSequence());
    if (seq_string.empty())
    {
      QMessageBox::warning(this, "Error", QString("You must enter a") + (peptide_input ? "peptide" : "RNA") + " sequence!");
      return;
    }
    
    AASequence aa_sequence;
    NASequence na_sequence;

    try
    {
      if (peptide_input)
      {
        aa_sequence = AASequence::fromString(seq_string);
      }
      else
      {
        na_sequence = NASequence::fromString(seq_string);
      }
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

    TheoreticalSpectrumGenerator pep_generator;
    NucleicAcidSpectrumGenerator na_generator;

    try
    {
      if (peptide_input)
      {
        pep_generator.setParameters(p);
        pep_generator.getSpectrum(spec_, aa_sequence, charge, charge);
      }
      else
      {
        na_generator.setParameters(p);
        na_generator.getSpectrum(spec_, na_sequence, charge, charge);
      }
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
    if (ui_.model_none->isChecked())
    {
      ui_.max_iso_label->setEnabled(false);
      ui_.max_iso_spinbox->setEnabled(false);
      ui_.max_iso_prob_label->setEnabled(false);
      ui_.max_iso_prob_spinbox->setEnabled(false);
    }
    else if (ui_.model_coarse->isChecked())
    {
      ui_.max_iso_label->setEnabled(true);
      ui_.max_iso_spinbox->setEnabled(true);
      ui_.max_iso_prob_label->setEnabled(false);
      ui_.max_iso_prob_spinbox->setEnabled(false);
    }
    else if (ui_.model_fine->isChecked())
    {
      ui_.max_iso_label->setEnabled(false);
      ui_.max_iso_spinbox->setEnabled(false);
      ui_.max_iso_prob_label->setEnabled(true);
      ui_.max_iso_prob_spinbox->setEnabled(true);
    }

  }

  void TheoreticalSpectrumGenerationDialog::seqTypeSwitch()
  {
    bool peptide_input = ui_.seq_type->currentText() == "Peptide";

    QListWidgetItem* a_b = ui_.list_widget->item(int(Checkbox::A_b_Ions));
    QListWidgetItem* d = ui_.list_widget->item(int(Checkbox::D_Ions));
    QListWidgetItem* w = ui_.list_widget->item(int(Checkbox::W_Ions));
    QListWidgetItem* losses = ui_.list_widget->item(int(Checkbox::Neutral_losses));
    QListWidgetItem* abundant_i = ui_.list_widget->item(int(Checkbox::Abundant_Immonium_Ions));

    if (peptide_input)
    {
      // enable isotopes
      ui_.isotope_model->setHidden(false);
      ui_.max_iso_label->setHidden(false);
      ui_.max_iso_spinbox->setHidden(false);
      ui_.max_iso_prob_label->setHidden(false);
      ui_.max_iso_prob_spinbox->setHidden(false);
      modelChanged();
      
      // ensable losses and immonium ions
      ui_.rel_loss_intensity->setHidden(false);
      ui_.rel_loss_label->setHidden(false);
      losses->setHidden(false);
      abundant_i->setHidden(false);

      // disable a-B-, D- and W-Ions

      ui_.a_b_intensity->setHidden(true);
      ui_.a_b_label->setHidden(true);
      a_b->setHidden(true);

      ui_.d_intensity->setHidden(true);
      ui_.d_label->setHidden(true);
      d->setHidden(true);

      ui_.w_intensity->setHidden(true);
      ui_.w_label->setHidden(true);
      w->setHidden(true);
    }
    else // rna input
    {
      // disable isotopes
      ui_.isotope_model->setHidden(true);
      ui_.max_iso_label->setHidden(true);
      ui_.max_iso_spinbox->setHidden(true);
      ui_.max_iso_prob_label->setHidden(true);
      ui_.max_iso_prob_spinbox->setHidden(true);

      // disable losses and immonium ions
      ui_.rel_loss_intensity->setHidden(true);
      ui_.rel_loss_label->setHidden(true);
      losses->setHidden(true);
      abundant_i->setHidden(true);

      // enable a-B-, D- and W-Ions

      ui_.a_b_intensity->setHidden(false);
      ui_.a_b_label->setHidden(false);
      a_b->setHidden(false);

      ui_.d_intensity->setHidden(false);
      ui_.d_label->setHidden(false);
      d->setHidden(false);

      ui_.w_intensity->setHidden(false);
      ui_.w_label->setHidden(false);
      w->setHidden(false);

      //this->resize(QSize(150, 300));
    }
  }

} // namespace

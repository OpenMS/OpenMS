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
#include <ui_TheoreticalSpectrumGenerationDialog.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/String.h>

// Qt includes
#include <QtWidgets/QMessageBox>
#include <qflags.h>

#include <array>

namespace OpenMS
{
  // for each check box (index) get corresponding parameter with description
  // Order is important here!
  // To access the right entry for each check box
  // use int(TheoreticalSpectrumGenerationDialog::CheckBox).
  const std::vector<std::pair<String, String>> check_box_to_param {
    {"add_a_ions", "Add peaks of a-ions to the spectrum"},
    {"add_a-B_ions", "Add peaks of a-B-ions to the spectrum (nucleotide sequences only)"},
    {"add_b_ions", "Add peaks of b-ions to the spectrum"},
    {"add_c_ions", "Add peaks of c-ions to the spectrum"},
    {"add_d_ions", "Add peaks of d-ions to the spectrum (nucleotide sequences only)"},
    {"add_w_ions", "Add peaks of w-ions to the spectrum (nucleotide sequences only)"},
    {"add_x_ions", "Add peaks of x-ions to the spectrum"},
    {"add_y_ions", "Add peaks of y-ions to the spectrum"},
    {"add_z_ions", "Add peaks of z-ions to the spectrum"},
    {"add_precursor_peaks", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes"},
    {"add_losses", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered (peptide sequences only)"},
    {"add_abundant_immonium_ions", "Add most abundant immonium ions (peptide sequences only)"}};

  TheoreticalSpectrumGenerationDialog::TheoreticalSpectrumGenerationDialog() : 
    ui_(new Ui::TheoreticalSpectrumGenerationDialogTemplate), check_box_to_intensity_{
                                                                  {&(ui_->a_intensity), "a"}, // A-Ion
                                                                  {&(ui_->a_b_intensity), "a-B"}, // a-B-Ion
                                                                  {&(ui_->b_intensity), "b"},   // B-Ion
                                                                  {&(ui_->c_intensity), "c"},   // C-Ion
                                                                  {&(ui_->d_intensity), "d"},   // D-Ion
                                                                  {&(ui_->w_intensity), "w"},   // W-Ion
                                                                  {&(ui_->x_intensity), "x"},   // X-Ion
                                                                  {&(ui_->y_intensity), "y"},   // Y-Ion
                                                                  {&(ui_->z_intensity), "z"},   // Z-Ion
                                                                  {nullptr, ""},            // Precursor
                                                                  {nullptr, ""},            // Neutral losses: ui_->rel_loss_intensity is a normal spin box
                                                                  {nullptr, ""}                // Abundant Immonium Ions
                                                                  },
      check_boxes_{CheckBox_(&(ui_->a_intensity), {State::ENABLED, State::ENABLED, State::HIDDEN}, {"add_a_ions", "Add peaks of a-ions to the spectrum"}, {"a_intensity", "Intensity of the a-ions"}),
                   CheckBox_(&(ui_->a_b_intensity), {State::HIDDEN, State::ENABLED, State::HIDDEN}, {"add_a-B_ions", "Add peaks of a-B-ions to the spectrum (nucleotide sequences only)"}, {"a-B_intensity", "Intensity of the a-B-ions"}),
                   CheckBox_(&(ui_->b_intensity), {State::PRECHECKED, State::ENABLED, State::HIDDEN}, {"add_b_ions", "Add peaks of b-ions to the spectrum"}, {"b_intensity", "Intensity of the b-ions"}),
                   CheckBox_(&(ui_->c_intensity), {State::ENABLED, State::ENABLED, State::HIDDEN}, {"add_c_ions", "Add peaks of c-ions to the spectrum"}, {"c_intensity", "Intensity of the c-ions"}),
                   CheckBox_(&(ui_->d_intensity), {State::HIDDEN, State::ENABLED, State::HIDDEN}, {"add_d_ions", "Add peaks of d-ions to the spectrum (nucleotide sequences only)"}, {"d_intensity", "Intensity of the d-ions"}),
                   CheckBox_(&(ui_->w_intensity), {State::HIDDEN, State::ENABLED, State::HIDDEN}, {"add_w_ions", "Add peaks of w-ions to the spectrum (nucleotide sequences only)"}, {"w_intensity", "Intensity of the w-ions"}),
                   CheckBox_(&(ui_->x_intensity), {State::ENABLED, State::ENABLED, State::HIDDEN}, {"add_x_ions", "Add peaks of x-ions to the spectrum"}, {"x_intensity", "Intensity of the x-ions"}),
                   CheckBox_(&(ui_->y_intensity), {State::PRECHECKED, State::ENABLED, State::HIDDEN}, {"add_y_ions", "Add peaks of y-ions to the spectrum"}, {"y_intensity", "Intensity of the y-ions"}),
                   CheckBox_(&(ui_->z_intensity), {State::ENABLED, State::ENABLED, State::HIDDEN}, {"add_z_ions", "Add peaks of z-ions to the spectrum"}, {"z_intensity", "Intensity of the z-ions"}),
                   CheckBox_(nullptr, {State::ENABLED, State::ENABLED, State::HIDDEN}, {"add_precursor_peaks", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes"}, {"", ""}),
                   // Neutral losses: ui_->rel_loss_intensity is a normal spin box
                   CheckBox_(nullptr, {State::ENABLED, State::HIDDEN, State::HIDDEN}, {"add_losses", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered (peptide sequences only)"}, {"", ""}),
                   CheckBox_(nullptr, {State::ENABLED, State::HIDDEN, State::HIDDEN}, {"add_abundant_immonium_ions", "Add most abundant immonium ions (peptide sequences only)"}, {"", ""})}
  {
    ui_->setupUi(this);
    
    // if dialog is accepted, try generating a spectrum, only close dialog on success
    connect(ui_->dialog_buttons, &QDialogButtonBox::accepted, this, &TheoreticalSpectrumGenerationDialog::calculateSpectrum_);

    // signals for changing isotope model interface
    connect(ui_->model_none_button, &QRadioButton::toggled, this, &TheoreticalSpectrumGenerationDialog::modelChanged_);
    connect(ui_->model_coarse_button, &QRadioButton::toggled, this, &TheoreticalSpectrumGenerationDialog::modelChanged_);
    connect(ui_->model_fine_button, &QRadioButton::toggled, this, &TheoreticalSpectrumGenerationDialog::modelChanged_);

    // for the list widget items are checked/unchecked if they are clicked on (disables clicking on the check box though ..)
    connect(ui_->ion_types, &QListWidget::itemClicked, this, &TheoreticalSpectrumGenerationDialog::listWidgetItemClicked_);

    // don't add any isotopes by default and update interface
    ui_->model_none_button->setChecked(true);
    modelChanged_();

    // signal for changing interface depending on sequence type
    connect(ui_->seq_type, &QComboBox::currentTextChanged, this, &TheoreticalSpectrumGenerationDialog::seqTypeSwitch_);

    // select peptide sequence by default and update interface
    ui_->seq_type->setCurrentText("Peptide");
    seqTypeSwitch_();

    // select b- and y-ions as residue types by default
    for (size_t i = 0; i < int(CheckBox::NUMBER_OF_CHECK_BOXES); ++i)
    {
      if (i == int(CheckBox::B_Ions) || i == int(CheckBox::Y_Ions))
      {
        ui_->ion_types->item(i)->setCheckState(Qt::Checked);
        continue;
      }
      ui_->ion_types->item(i)->setCheckState(Qt::Unchecked);
    }

    // automatic layout
    // disables manual resizing from the user
    layout()->setSizeConstraint(QLayout::SetFixedSize);
  }

  TheoreticalSpectrumGenerationDialog::~TheoreticalSpectrumGenerationDialog()
  {
    delete ui_;
  }

  const String TheoreticalSpectrumGenerationDialog::getSequence() const
  {
    return ui_->seq_input->text();
  }

  Param TheoreticalSpectrumGenerationDialog::getParam_() const
  {
    Param p;

    if (seq_type_ != "Metabolite") // no ions for metabolite input
    {
      // add check boxes to parameters, i.e. ion types
      for (size_t i = 0; i < int(CheckBox::NUMBER_OF_CHECK_BOXES); ++i)
      {
        // for peptide input skip rna specific ions
        if (seq_type_ == "Peptide" && (check_boxes_.at(i).state.at(0) == State::HIDDEN))
          continue;

        // for rna input skip peptide specific ions
        if (seq_type_ == "RNA" && (check_boxes_.at(i).state.at(1) == State::HIDDEN))
          continue;

        // set ion itself
        bool status = (ui_->ion_types->item(i)->checkState() == Qt::Checked);
        p.setValue(check_boxes_.at(i).param_this.first, status ? "true" : "false", check_boxes_.at(i).param_this.second);

        // set intensity of ion
        if (status)
        {
          QDoubleSpinBox** spin_ptr = check_boxes_.at(i).ptr_to_spin_box;
          if (spin_ptr == nullptr)
            continue;

          p.setValue(check_boxes_.at(i).param_spin.first, (*spin_ptr)->value(), check_boxes_.at(i).param_spin.second);
        }
      }
    }

    // charge
    p.setValue("charge", ui_->charge_spinbox->value());

    if (!(seq_type_ == "RNA")) // skip isotope pattern for RNA input
    {
      // isotopes
      if (!ui_->model_none_button->isChecked()) // add isotopes if any other model than 'None' is chosen
      {
        bool coarse_model = ui_->model_coarse_button->isChecked();
        String model = coarse_model ? "coarse" : "fine";
        p.setValue("isotope_model", model,
                   "Model to use for isotopic peaks ('none' means no isotopic peaks are added, 'coarse' adds isotopic peaks in unit mass distance, 'fine' uses the hyperfine isotopic generator to add "
                   "accurate isotopic peaks. Note that adding isotopic peaks is very slow.");

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
        p.setValue("isotope_model", "none",
                   "Model to use for isotopic peaks ('none' means no isotopic peaks are added, 'coarse' adds isotopic peaks in unit mass distance, 'fine' uses the hyperfine isotopic generator to add "
                   "accurate isotopic peaks. Note that adding isotopic peaks is very slow.");
      }

      if (seq_type_ == "Peptide")
      {
        // loss intensity
        double rel_loss_int = (double)(ui_->rel_loss_intensity->value()) / 100.0;
        p.setValue("relative_loss_intensity", rel_loss_int, "Intensity of loss ions, in relation to the intact ion intensity");
      }
    }

    return p;
  }

  const MSSpectrum& TheoreticalSpectrumGenerationDialog::getSpectrum() const
  {
    return spec_;
  }

  void TheoreticalSpectrumGenerationDialog::calculateSpectrum_()
  {
    if (!spec_.empty()) spec_.clear(true);

    String seq_string(this->getSequence());
    if (seq_string.empty())
    {
      QMessageBox::warning(this, "Error", QString("You must enter a ") + seq_type_ + " sequence!");
      return;
    }
    
    AASequence aa_sequence;
    NASequence na_sequence;
    EmpiricalFormula ef;

    try
    {
      if (seq_type_ == "Peptide")
      {
        aa_sequence = AASequence::fromString(seq_string);
      }
      else if (seq_type_ == "RNA")
      {
        na_sequence = NASequence::fromString(seq_string);
      }
      else
      {
        ef = EmpiricalFormula::fromString(seq_string);
      }
    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::warning(this, "Error", QString("Spectrum generation failed! (") + e.what() + ")");
      return;
    }

    Param p = this->getParam_();
    Int charge = p.getValue("charge");
    p.remove("charge"); // "charge" isn't a parameter of TheoreticalSpectrumGenerator

    p.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");

    TheoreticalSpectrumGenerator pep_generator;
    NucleicAcidSpectrumGenerator na_generator;

    try
    {
      if (seq_type_ == "Peptide")
      {
        pep_generator.setParameters(p);
        pep_generator.getSpectrum(spec_, aa_sequence, charge, charge);
      }
      else if (seq_type_ == "RNA")
      {
        na_generator.setParameters(p);
        na_generator.getSpectrum(spec_, na_sequence, charge, charge);
      }
      else // metabolite
      {
        IsotopeDistribution dist;
        if (p.getValue("isotope_model") == "coarse")
        {
          dist = ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(Int(p.getValue("max_isotope"))));
        }
        else if (p.getValue("isotope_model") == "fine")
        {
          dist = ef.getIsotopeDistribution(FineIsotopePatternGenerator(double(p.getValue("max_isotope_probability"))));
        }
        else
        {
          QMessageBox::warning(this, "Error", QString("Isotope model 'None' is not supported for metabolite input!"));
          return;
        }
        for (const auto& it : dist)
        {
          // currently no meta data is written in this setting
          spec_.emplace_back(it.getMZ() / charge, it.getIntensity());
        }
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

  void TheoreticalSpectrumGenerationDialog::modelChanged_()
  {    
    if (ui_->model_none_button->isChecked())
    {
      ui_->max_iso_label->setEnabled(false);
      ui_->max_iso_spinbox->setEnabled(false);
      ui_->max_iso_prob_label->setEnabled(false);
      ui_->max_iso_prob_spinbox->setEnabled(false);
    }
    else if (ui_->model_coarse_button->isChecked())
    {
      ui_->max_iso_label->setEnabled(true);
      ui_->max_iso_spinbox->setEnabled(true);
      ui_->max_iso_prob_label->setEnabled(false);
      ui_->max_iso_prob_spinbox->setEnabled(false);
    }
    else if (ui_->model_fine_button->isChecked())
    {
      ui_->max_iso_label->setEnabled(false);
      ui_->max_iso_spinbox->setEnabled(false);
      ui_->max_iso_prob_label->setEnabled(true);
      ui_->max_iso_prob_spinbox->setEnabled(true);
    }

  }

  void TheoreticalSpectrumGenerationDialog::seqTypeSwitch_()
  {
    seq_type_ = ui_->seq_type->currentText(); // save current seq type setting in member

    if (seq_type_ == "Peptide")
    {
      // enable ion types and intensities
      ui_->ion_types->setHidden(false);
      ui_->ion_types_label->setHidden(false);
      ui_->intensities->setHidden(false);

      // enable isotopes
      ui_->isotope_model->setHidden(false);
      ui_->max_iso_label->setHidden(false);
      ui_->max_iso_spinbox->setHidden(false);
      ui_->max_iso_prob_label->setHidden(false);
      ui_->max_iso_prob_spinbox->setHidden(false);
      modelChanged_();

      // enable isotope model 'none'
      ui_->model_none_button->setEnabled(true);
      
      // ensable losses and immonium ions
      ui_->rel_loss_intensity->setHidden(false);
      ui_->rel_loss_label->setHidden(false);
      losses->setHidden(false);
      abundant_i->setHidden(false);

      // disable a-B-, D- and W-Ions

      ui_->a_b_intensity->setHidden(true);
      ui_->a_b_label->setHidden(true);
      a_b->setHidden(true);

      ui_->d_intensity->setHidden(true);
      ui_->d_label->setHidden(true);
      d->setHidden(true);

      ui_->w_intensity->setHidden(true);
      ui_->w_label->setHidden(true);
      w->setHidden(true);
    }
    else
    {
      if (seq_type_ == "RNA") // rna input
      {
        // enable ion types and intensities
        ui_->ion_types->setHidden(false);
        ui_->ion_types_label->setHidden(false);
        ui_->intensities->setHidden(false);

        // disable isotopes
        ui_->isotope_model->setHidden(true);
        ui_->max_iso_label->setHidden(true);
        ui_->max_iso_spinbox->setHidden(true);
        ui_->max_iso_prob_label->setHidden(true);
        ui_->max_iso_prob_spinbox->setHidden(true);

        // disable losses and immonium ions
        ui_->rel_loss_intensity->setHidden(true);
        ui_->rel_loss_label->setHidden(true);
        losses->setHidden(true);
        abundant_i->setHidden(true);

        // enable a-B-, D- and W-Ions

        ui_->a_b_intensity->setHidden(false);
        ui_->a_b_label->setHidden(false);
        a_b->setHidden(false);

        ui_->d_intensity->setHidden(false);
        ui_->d_label->setHidden(false);
        d->setHidden(false);

        ui_->w_intensity->setHidden(false);
        ui_->w_label->setHidden(false);
        w->setHidden(false);
      }
      else // metabolite input
      {
        // enable isotopes
        ui_->isotope_model->setHidden(false);
        ui_->max_iso_label->setHidden(false);
        ui_->max_iso_spinbox->setHidden(false);
        ui_->max_iso_prob_label->setHidden(false);
        ui_->max_iso_prob_spinbox->setHidden(false);
        modelChanged_();
        
        // disable isotope model 'none'
        ui_->model_none_button->setEnabled(false);
        if (ui_->model_none_button->isChecked())
        {
          ui_->model_fine_button->setChecked(Qt::CheckState::Checked);
        }

        // disable ion types and intensities
        ui_->ion_types->setHidden(true);
        ui_->ion_types_label->setHidden(true);
        ui_->intensities->setHidden(true);
      }
    }
  }

  void TheoreticalSpectrumGenerationDialog::listWidgetItemClicked_(QListWidgetItem* item)
  {
    if (item->checkState() == Qt::CheckState::Checked)
    {
      item->setCheckState(Qt::CheckState::Unchecked);
      return;
    }
    item->setCheckState(Qt::CheckState::Checked);
    return;
  }
  
  TheoreticalSpectrumGenerationDialog::CheckBox_::CheckBox_(QDoubleSpinBox** p, std::array<State, 3> s, std::pair<String, String> p_t, std::pair<String, String> p_s)
  {
    ptr_to_spin_box = p;
    state = s;
    param_this = p_t;
    param_spin = p_s;
  }
} // namespace OpenMS

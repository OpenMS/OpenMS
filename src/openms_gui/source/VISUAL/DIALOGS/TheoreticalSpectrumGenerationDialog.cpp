// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <utility>

namespace OpenMS
{
  TheoreticalSpectrumGenerationDialog::CheckBox::CheckBox(QDoubleSpinBox** sb, QLabel** l, std::array<CheckBoxState, 3> s, std::pair<String, String> p_t, std::pair<String, String> p_s) :
      ptr_to_spin_box(sb), ptr_to_spin_label(l), state(s), param_this(std::move(p_t)), param_spin(std::move(p_s))
  {}

  TheoreticalSpectrumGenerationDialog::TheoreticalSpectrumGenerationDialog() : 
    ui_(new Ui::TheoreticalSpectrumGenerationDialogTemplate),

    // Order has to be the same as in the UI!
    check_boxes_ {
        CheckBox(&(ui_->a_intensity), &(ui_->a_label), {CheckBoxState::ENABLED, CheckBoxState::ENABLED, CheckBoxState::HIDDEN}, {"add_a_ions", "Add peaks of a-ions to the spectrum"},
                  {"a_intensity", "Intensity of the a-ions"}),
        CheckBox(&(ui_->a_b_intensity), &(ui_->a_b_label), {CheckBoxState::HIDDEN, CheckBoxState::ENABLED, CheckBoxState::HIDDEN}, {"add_a-B_ions", "Add peaks of a-B-ions to the spectrum (nucleotide sequences only)"},
                  {"a-B_intensity", "Intensity of the a-B-ions"}),
        CheckBox(&(ui_->b_intensity), &(ui_->b_label), {CheckBoxState::PRECHECKED, CheckBoxState::PRECHECKED, CheckBoxState::HIDDEN}, {"add_b_ions", "Add peaks of b-ions to the spectrum"},
                  {"b_intensity", "Intensity of the b-ions"}),
        CheckBox(&(ui_->c_intensity), &(ui_->c_label), {CheckBoxState::ENABLED, CheckBoxState::ENABLED, CheckBoxState::HIDDEN}, {"add_c_ions", "Add peaks of c-ions to the spectrum"},
                  {"c_intensity", "Intensity of the c-ions"}),
        CheckBox(&(ui_->d_intensity), &(ui_->d_label), {CheckBoxState::HIDDEN, CheckBoxState::ENABLED, CheckBoxState::HIDDEN}, {"add_d_ions", "Add peaks of d-ions to the spectrum (nucleotide sequences only)"},
                  {"d_intensity", "Intensity of the d-ions"}),
        CheckBox(&(ui_->w_intensity), &(ui_->w_label), {CheckBoxState::HIDDEN, CheckBoxState::ENABLED, CheckBoxState::HIDDEN}, {"add_w_ions", "Add peaks of w-ions to the spectrum (nucleotide sequences only)"},
                  {"w_intensity", "Intensity of the w-ions"}),
        CheckBox(&(ui_->x_intensity), &(ui_->x_label), {CheckBoxState::ENABLED, CheckBoxState::ENABLED, CheckBoxState::HIDDEN}, {"add_x_ions", "Add peaks of x-ions to the spectrum"},
                  {"x_intensity", "Intensity of the x-ions"}),
        CheckBox(&(ui_->y_intensity), &(ui_->y_label), {CheckBoxState::PRECHECKED, CheckBoxState::PRECHECKED, CheckBoxState::HIDDEN}, {"add_y_ions", "Add peaks of y-ions to the spectrum"},
                  {"y_intensity", "Intensity of the y-ions"}),
        CheckBox(&(ui_->z_intensity), &(ui_->z_label), {CheckBoxState::ENABLED, CheckBoxState::ENABLED, CheckBoxState::HIDDEN}, {"add_z_ions", "Add peaks of z-ions to the spectrum"},
                  {"z_intensity", "Intensity of the z-ions"}),
        CheckBox(nullptr, nullptr, {CheckBoxState::ENABLED, CheckBoxState::ENABLED, CheckBoxState::HIDDEN},
                              {"add_precursor_peaks", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes"}, {"", ""}),
        // Neutral losses: ui_->rel_loss_intensity is a normal spin box and has to be checked manually
        CheckBox(nullptr, nullptr, {CheckBoxState::ENABLED, CheckBoxState::HIDDEN, CheckBoxState::HIDDEN},
                  {"add_losses", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered (peptide sequences only)"}, {"", ""}),
        CheckBox(nullptr, nullptr, {CheckBoxState::ENABLED, CheckBoxState::HIDDEN, CheckBoxState::HIDDEN},
                  {"add_abundant_immonium_ions", "Add most abundant immonium ions (peptide sequences only)"}, {"", ""})}
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
    modelChanged_(); //because setting the check state of the button doesn't call 'toggled'

    // signal for changing interface depending on sequence type
    connect(ui_->seq_type, &QComboBox::currentTextChanged, this, &TheoreticalSpectrumGenerationDialog::seqTypeSwitch_);

    // To set the interface and members
    seqTypeSwitch_();

    
    for (size_t i = 0; i < check_boxes_.size(); ++i)
    {
      if (check_boxes_.at(i).state.at(0) == CheckBoxState::PRECHECKED) // 'state.at(0)' because seq type was just set to "Peptide"
      {
        ui_->ion_types->item(i)->setCheckState(Qt::Checked);
      }
      else
      {
        ui_->ion_types->item(i)->setCheckState(Qt::Unchecked);
      }
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

    if (seq_type_ != SequenceType::METABOLITE) // no ions for metabolite input
    {
      // add check boxes to parameters, i.e. ion types
      for (size_t i = 0; i < check_boxes_.size(); ++i)
      {
        // for peptide input skip rna specific ions
        if (seq_type_ == SequenceType::PEPTIDE && (check_boxes_.at(i).state.at(0) == CheckBoxState::HIDDEN))
          continue;

        // for rna input skip peptide specific ions
        if (seq_type_ == SequenceType::RNA && (check_boxes_.at(i).state.at(1) == CheckBoxState::HIDDEN))
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

    if (!(seq_type_ == SequenceType::RNA)) // skip isotope pattern for RNA input
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

      if (seq_type_ == SequenceType::PEPTIDE)
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
      const std::array<String, 3> types{"Peptide", "RNA", "Metabolite"};
      QMessageBox::warning(this, "Error", QString("You must enter a ") + QString::fromStdString(types.at(int(seq_type_))) + " sequence!");
      return;
    }
    
    AASequence aa_sequence;
    NASequence na_sequence;
    EmpiricalFormula ef;

    try
    {
      if (seq_type_ == SequenceType::PEPTIDE)
      {
        aa_sequence = AASequence::fromString(seq_string);
      }
      else if (seq_type_ == SequenceType::RNA)
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
      if (seq_type_ == SequenceType::PEPTIDE)
      {
        p.setValue("add_first_prefix_ion", "true"); // do not skip b1 ion
        pep_generator.setParameters(p);
        pep_generator.getSpectrum(spec_, aa_sequence, charge, charge);
      }
      else if (seq_type_ == SequenceType::RNA)
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
    // save current sequence type setting in member
    String tmp = ui_->seq_type->currentText();
    if (tmp == "Peptide")
    {
      seq_type_ = SequenceType::PEPTIDE;
    }
    else if (tmp == "RNA")
    {
      seq_type_ = SequenceType::RNA;
    }
    else if (tmp == "Metabolite")
    {
      seq_type_ = SequenceType::METABOLITE;
    }
    else
    {
      // this will be reached if the entries of the sequence type combo box are changed
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Setting for the sequence type unkown! Was: " + tmp);
    }

    if (seq_type_ == SequenceType::PEPTIDE)
    {
      // change sequence label
      ui_->enter_seq_label->setText("Enter sequence: ");

      // enable ion types and intensities
      ui_->ion_types->setHidden(false);
      ui_->ion_types_label->setHidden(false);
      ui_->intensities->setHidden(false);
      updateIonTypes_();

      // enable isotopes
      ui_->isotope_model->setHidden(false);
      modelChanged_();

      // enable isotope model 'none'
      ui_->model_none_button->setEnabled(true);
    }
    else
    {
      if (seq_type_ == SequenceType::RNA) // rna input
      {
        // change sequence label
        ui_->enter_seq_label->setText("Enter sequence: ");

        // enable ion types and intensities
        ui_->ion_types->setHidden(false);
        ui_->ion_types_label->setHidden(false);
        ui_->intensities->setHidden(false);
        updateIonTypes_();

        // disable isotopes
        ui_->isotope_model->setHidden(true);
      }
      else // metabolite input
      {
        // change sequence label
        ui_->enter_seq_label->setText("Enter empirical formula (e.g. C6H12O6): ");

        // enable isotopes
        ui_->isotope_model->setHidden(false);
        modelChanged_();
        
        // disable isotope model 'none'
        ui_->model_none_button->setEnabled(false);
        if (ui_->model_none_button->isChecked())
        {
          ui_->model_fine_button->setChecked(true);
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
  
  void TheoreticalSpectrumGenerationDialog::updateIonTypes_()
  {
    int input_type;
    if (seq_type_ == SequenceType::PEPTIDE)
      input_type = 0;
    else if (seq_type_ == SequenceType::RNA)
      input_type = 1;
    else // Metabolite
      input_type = 2;

    for (size_t i = 0; i < check_boxes_.size(); ++i)
    {
      const CheckBox* curr_box = &check_boxes_.at(i);

      bool hidden(curr_box->state.at(input_type) == CheckBoxState::HIDDEN);

      // activate check box
      ui_->ion_types->item(i)->setHidden(hidden);
      
      // activte intensity with label
      QDoubleSpinBox** spin_ptr = curr_box->ptr_to_spin_box;
      if (spin_ptr == nullptr)
      {
        // manually check for neutral losses
        if (curr_box->param_this.first == "add_losses")
        {
          ui_->rel_loss_intensity->setHidden(hidden);
          ui_->rel_loss_label->setHidden(hidden);
        }
        continue;
      }
      (*spin_ptr)->setHidden(hidden);

      QLabel** label_ptr = curr_box->ptr_to_spin_label;
      if (label_ptr == nullptr)
        continue;
      (*label_ptr)->setHidden(hidden);
    }
  }
} // namespace OpenMS

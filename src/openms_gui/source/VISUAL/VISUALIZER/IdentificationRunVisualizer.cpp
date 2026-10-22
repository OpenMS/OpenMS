// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/IdentificationRunVisualizer.h>

#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

//QT
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QTextEdit>

namespace OpenMS
{

  IdentificationRunVisualizer::IdentificationRunVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<IdentificationData::RunSettings>()
  {
    addLabel_("Modify the settings of the identification run.");
    addSeparator_();
    addLineEdit_(identifier_, "Identifier");
    addLineEdit_(kind_, "Molecule kind");
    addTextEdit_(databases_, "Databases");
    addTextEdit_(scores_, "Score definitions<br>(of all runs)");
    for (auto* field : {identifier_, kind_}) field->setReadOnly(true);
    for (auto* field : {databases_, scores_}) field->setReadOnly(true);
    addSeparator_();
    addLineEdit_(software_, "Search engine");
    addLineEdit_(software_version_, "Search engine version");
    addLineEdit_(date_, "Date of search");
    addSeparator_();
    addLabel_("Search Parameters:");
    addLineEdit_(charges_, "Allowed charges");
    addIntLineEdit_(missed_cleavages_, "Missed Cleavages");
    addDoubleLineEdit_(fragment_tolerance_, "Fragment ion mass tolerance");
    addBooleanComboBox_(fragment_tolerance_ppm_, "Fragment ion mass tolerance in ppm");
    addDoubleLineEdit_(precursor_tolerance_, "Precursor ion mass tolerance");
    addBooleanComboBox_(precursor_tolerance_ppm_, "Precursor ion mass tolerance in ppm");
    addComboBox_(mass_type_, "Mass type");
    addLineEdit_(enzyme_, "Digestion enzyme");
    // (the object names identify the fields, e.g. in tests)
    software_->setObjectName("software");
    software_version_->setObjectName("software_version");
    date_->setObjectName("date");
    charges_->setObjectName("charges");

    finishAdding_();
  }

  void IdentificationRunVisualizer::setRun(const IdentificationData::Run& run)
  {
    identifier_->setText(toQString(run.getIdentifier()));
    switch (run.getMoleculeKind())
    {
      case IdentificationData::MoleculeKind::PEPTIDE: kind_->setText("peptide"); break;
      case IdentificationData::MoleculeKind::OLIGONUCLEOTIDE: kind_->setText("oligonucleotide"); break;
      case IdentificationData::MoleculeKind::COMPOUND: kind_->setText("compound"); break;
    }
    QStringList databases;
    for (const auto& database : run.getDatabases())
    {
      QString text = toQString(database.path);
      if (!database.version.empty()) text += " (version " + toQString(database.version) + ")";
      if (!database.taxonomy.empty()) text += ", taxonomy " + toQString(database.taxonomy);
      databases << text;
    }
    databases_->setText(databases.join("\n"));
    QStringList scores;
    const auto primary = run.getPrimaryScore();
    for (Size i = 0; i < run.getScoreDefinitions().size(); ++i)
    {
      const auto& definition = run.getScoreDefinitions()[i];
      QString text = toQString(definition.name) + (definition.higher_better ? " (higher is better)" : " (lower is better)");
      if (primary && primary->value == i) text += ", primary";
      scores << text;
    }
    scores_->setText(scores.join("\n"));
  }

  void IdentificationRunVisualizer::update_()
  {
    software_->setText(toQString(temp_.software));
    software_version_->setText(toQString(temp_.software_version));
    date_->setText(temp_.date.isValid() ? toQString(temp_.date.get()) : QString());
    const auto& search = temp_.search;
    charges_->setText(toQString(search.charges));
    missed_cleavages_->setText(QString::number(search.missed_cleavages));
    fragment_tolerance_->setText(QString::number(search.fragment_mass_tolerance));
    fragment_tolerance_ppm_->setCurrentIndex(search.fragment_mass_tolerance_ppm);
    precursor_tolerance_->setText(QString::number(search.precursor_mass_tolerance));
    precursor_tolerance_ppm_->setCurrentIndex(search.precursor_mass_tolerance_ppm);
    enzyme_->setText(toQString(search.digestion_enzyme.getName()));
    const int n_mass_types = static_cast<int>(ProteinIdentification::PeakMassType::SIZE_OF_PEAKMASSTYPE);
    if (!isEditable())
    {
      fillComboBox_(mass_type_, &ProteinIdentification::NamesOfPeakMassType[static_cast<size_t>(search.mass_type)], 1);
    }
    else
    {
      fillComboBox_(mass_type_, ProteinIdentification::NamesOfPeakMassType, n_mass_types);
      mass_type_->setCurrentIndex(static_cast<int>(search.mass_type));
    }
  }

  void IdentificationRunVisualizer::store()
  {
    // a field is written only if its text changed, so values that the text shows rounded stay exact
    const auto& search = temp_.search;
    const auto changed = [](const QLineEdit* field, const QString& loaded) { return field->text() != loaded; };
    if (changed(software_, toQString(temp_.software))) ptr_->software = fromQString(software_->text());
    if (changed(software_version_, toQString(temp_.software_version))) ptr_->software_version = fromQString(software_version_->text());
    if (changed(date_, temp_.date.isValid() ? toQString(temp_.date.get()) : QString()))
    {
      const auto date_text = fromQString(date_->text());
      if (date_text.empty())
      {
        ptr_->date = DateTime();
      }
      else
      {
        try
        {
          DateTime date;
          date.set(date_text);
          ptr_->date = date;
        }
        catch (std::exception&)
        {
          emit sendStatus("Format of the date of the identification run is not correct.");
        }
      }
    }
    if (changed(charges_, toQString(search.charges))) ptr_->search.charges = fromQString(charges_->text());
    if (changed(missed_cleavages_, QString::number(search.missed_cleavages))) ptr_->search.missed_cleavages = missed_cleavages_->text().toInt();
    if (changed(fragment_tolerance_, QString::number(search.fragment_mass_tolerance)))
    {
      ptr_->search.fragment_mass_tolerance = fragment_tolerance_->text().toDouble();
    }
    ptr_->search.fragment_mass_tolerance_ppm = fragment_tolerance_ppm_->currentIndex();
    if (changed(precursor_tolerance_, QString::number(search.precursor_mass_tolerance)))
    {
      ptr_->search.precursor_mass_tolerance = precursor_tolerance_->text().toDouble();
    }
    ptr_->search.precursor_mass_tolerance_ppm = precursor_tolerance_ppm_->currentIndex();
    if (isEditable()) ptr_->search.mass_type = static_cast<ProteinIdentification::PeakMassType>(mass_type_->currentIndex());
    const auto enzyme = fromQString(enzyme_->text());
    if (enzyme != search.digestion_enzyme.getName())
    {
      if (ProteaseDB::getInstance()->hasEnzyme(enzyme))
      {
        ptr_->search.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme));
      }
      else
      {
        emit sendStatus("Unknown digestion enzyme of the identification run: " + enzyme);
      }
    }

    temp_ = (*ptr_);
  }

  void IdentificationRunVisualizer::undo_()
  {
    update_();
  }

}

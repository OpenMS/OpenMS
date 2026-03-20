// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/VISUAL/VISUALIZER/ProteinIdentificationVisualizer.h>
#include <OpenMS/VISUAL/MISC/QtHelpers.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
//QT
#include <QtWidgets/QLineEdit>
#include <QValidator>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QComboBox>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{
  ProteinIdentificationVisualizer::ProteinIdentificationVisualizer(bool editable, QWidget * parent, MetaDataBrowser * caller) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<ProteinIdentification>()
  {
    pidv_caller_ = caller;

    addLineEdit_(identifier_, "Identifier<br>(of corresponding PeptideIdentifications)");

    addSeparator_();
    addLineEdit_(engine_, "Search engine");
    addLineEdit_(engine_version_, "Search engine version");
    addLineEdit_(identification_date_, "Date of search");
    addLineEdit_(score_type_, "Score type");
    addBooleanComboBox_(higher_better_, "Higher score is better");
    addDoubleLineEdit_(identification_threshold_, "Protein significance threshold");

    addSeparator_();
    addLabel_("Search Parameters:");
    addLineEdit_(db_, "Database name");
    addLineEdit_(db_version_, "Database version");
    addLineEdit_(taxonomy_, "Taxonomy restriction");
    addLineEdit_(charges_, "Allowed charges");
    addIntLineEdit_(missed_cleavages_, "Missed Cleavages");
    addDoubleLineEdit_(peak_tolerance_, "Fragment ion mass tolerance");
    addDoubleLineEdit_(precursor_tolerance_, "Precursor ion mass tolerance");
    addComboBox_(mass_type_, "Mass type");
    addLineEdit_(enzyme_, "Digestion enzyme");

    addSeparator_();
    addLabel_("Show protein hits with score equal or better than a threshold.");
    QPushButton * button;
    addLineEditButton_("Score threshold", filter_threshold_, button, "Filter");
    connect(button, SIGNAL(clicked()), this, SLOT(updateTree_()));

    finishAdding_();
  }

  void ProteinIdentificationVisualizer::load(ProteinIdentification & s, int tree_item_id)
  {
    ptr_ = &s;
    temp_ = s;

    // id of the item in the tree
    tree_id_ = tree_item_id;

    identification_date_->setText(toQString(temp_.getDateTime().get()));
    identification_threshold_->setText(QString::number(temp_.getSignificanceThreshold()));
    identifier_->setText(toQString(temp_.getIdentifier()));
    engine_->setText(toQString(temp_.getSearchEngine()));
    engine_version_->setText(toQString(temp_.getSearchEngineVersion()));
    score_type_->setText(toQString(temp_.getScoreType()));
    higher_better_->setCurrentIndex(temp_.isHigherScoreBetter());

    db_->setText(toQString(temp_.getSearchParameters().db));
    db_version_->setText(toQString(temp_.getSearchParameters().db_version));
    taxonomy_->setText(toQString(temp_.getSearchParameters().taxonomy));
    charges_->setText(toQString(temp_.getSearchParameters().charges));
    missed_cleavages_->setText(QString::number(temp_.getSearchParameters().missed_cleavages));
    peak_tolerance_->setText(QString::number(temp_.getSearchParameters().fragment_mass_tolerance));
    precursor_tolerance_->setText(QString::number(temp_.getSearchParameters().precursor_mass_tolerance));
    enzyme_->setText(toQString(temp_.getSearchParameters().digestion_enzyme.getName()));

    if (!isEditable())
    {
      fillComboBox_(mass_type_, &ProteinIdentification::NamesOfPeakMassType[static_cast<size_t>(temp_.getSearchParameters().mass_type)], 1);
    }
    else
    {
      fillComboBox_(mass_type_, ProteinIdentification::NamesOfPeakMassType, static_cast<int>(ProteinIdentification::PeakMassType::SIZE_OF_PEAKMASSTYPE));
      mass_type_->setCurrentIndex(static_cast<int>(temp_.getSearchParameters().mass_type));
    }
  }

  void ProteinIdentificationVisualizer::updateTree_()
  {
    if (filter_threshold_->text() != "")
    {
      pidv_caller_->filterHits_(filter_threshold_->text().toDouble(), temp_.isHigherScoreBetter(), tree_id_);
    }
    else
    {
      pidv_caller_->showAllHits_(tree_id_);
    }
  }

  void ProteinIdentificationVisualizer::store()
  {
    ptr_->setSearchEngine(engine_->text().toStdString());
    ptr_->setSearchEngineVersion(engine_version_->text().toStdString());
    ptr_->setIdentifier(identifier_->text().toStdString());
    ptr_->setSignificanceThreshold(identification_threshold_->text().toFloat());
    ptr_->setScoreType(score_type_->text().toStdString());
    ptr_->setHigherScoreBetter(higher_better_->currentIndex());
    //date
    DateTime date;
    try
    {
      date.set(identification_date_->text().toStdString());
      ptr_->setDateTime(date);
    }
    catch (exception & /*e*/)
    {
      if (date.isNull())
      {
        std::string status = "Format of date in PROTEINIDENTIFICATION is not correct.";
        emit sendStatus(status);
      }
    }

    //search parameters
    ProteinIdentification::SearchParameters tmp = ptr_->getSearchParameters();
    tmp.db = db_->text().toStdString();
    tmp.db_version = db_version_->text().toStdString();
    tmp.taxonomy = taxonomy_->text().toStdString();
    tmp.charges = charges_->text().toStdString();
    tmp.missed_cleavages = missed_cleavages_->text().toInt();
    tmp.fragment_mass_tolerance = peak_tolerance_->text().toFloat();
    tmp.precursor_mass_tolerance = precursor_tolerance_->text().toFloat();
    tmp.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme_->text().toStdString()));
    tmp.mass_type = (ProteinIdentification::PeakMassType)(mass_type_->currentIndex());
    ptr_->setSearchParameters(tmp);

    temp_ = (*ptr_);
  }

  void ProteinIdentificationVisualizer::undo_()
  {
    load(*ptr_, tree_id_);
  }

}

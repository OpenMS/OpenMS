// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/VISUAL/VISUALIZER/ProteinIdentificationVisualizer.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
//QT
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>
#include <QtGui/QPushButton>
#include <QtGui/QComboBox>

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

    identification_date_->setText(temp_.getDateTime().get().toQString());
    identification_threshold_->setText(QString::number(temp_.getSignificanceThreshold()));
    identifier_->setText(temp_.getIdentifier().toQString());
    engine_->setText(temp_.getSearchEngine().toQString());
    engine_version_->setText(temp_.getSearchEngineVersion().toQString());
    score_type_->setText(temp_.getScoreType().toQString());
    higher_better_->setCurrentIndex(temp_.isHigherScoreBetter());

    db_->setText(temp_.getSearchParameters().db.toQString());
    db_version_->setText(temp_.getSearchParameters().db_version.toQString());
    taxonomy_->setText(temp_.getSearchParameters().taxonomy.toQString());
    charges_->setText(temp_.getSearchParameters().charges.toQString());
    missed_cleavages_->setText(QString::number(temp_.getSearchParameters().missed_cleavages));
    peak_tolerance_->setText(QString::number(temp_.getSearchParameters().fragment_mass_tolerance));
    precursor_tolerance_->setText(QString::number(temp_.getSearchParameters().precursor_mass_tolerance));
    enzyme_->setText(temp_.getSearchParameters().digestion_enzyme.getName().toQString());

    if (!isEditable())
    {
      fillComboBox_(mass_type_, &ProteinIdentification::NamesOfPeakMassType[temp_.getSearchParameters().mass_type], 1);
    }
    else
    {
      fillComboBox_(mass_type_, ProteinIdentification::NamesOfPeakMassType, ProteinIdentification::SIZE_OF_PEAKMASSTYPE);
      mass_type_->setCurrentIndex(temp_.getSearchParameters().mass_type);
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
    ptr_->setSearchEngine(engine_->text());
    ptr_->setSearchEngineVersion(engine_version_->text());
    ptr_->setIdentifier(identifier_->text());
    ptr_->setSignificanceThreshold(identification_threshold_->text().toFloat());
    ptr_->setScoreType(score_type_->text());
    ptr_->setHigherScoreBetter(higher_better_->currentIndex());
    //date
    DateTime date;
    try
    {
      date.set(identification_date_->text());
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
    tmp.db = db_->text();
    tmp.db_version = db_version_->text();
    tmp.taxonomy = taxonomy_->text();
    tmp.charges = charges_->text();
    tmp.missed_cleavages = missed_cleavages_->text().toInt();
    tmp.fragment_mass_tolerance = peak_tolerance_->text().toFloat();
    tmp.precursor_mass_tolerance = precursor_tolerance_->text().toFloat();
    tmp.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme_->text()));
    tmp.mass_type = (ProteinIdentification::PeakMassType)(mass_type_->currentIndex());
    ptr_->setSearchParameters(tmp);

    temp_ = (*ptr_);
  }

  void ProteinIdentificationVisualizer::undo_()
  {
    load(*ptr_, tree_id_);
  }

}

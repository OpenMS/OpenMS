// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/SwathLibraryStats.h>
#include <ui_SwathLibraryStats.h>

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

using namespace std;

namespace OpenMS
{
   
  SwathLibraryStats::SwathLibraryStats(QWidget* parent) :
      QWidget(parent),
      ui_(new Ui::SwathLibraryStats)
  {
    ui_->setupUi(this);
    ui_->table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeMode::ResizeToContents);
  }

  SwathLibraryStats::~SwathLibraryStats()
  {
    delete ui_;
  }

  void SwathLibraryStats::update(const TargetedExperiment::SummaryStatistics& stats)
  {
    auto getItem = [](const QString& text) {
      auto item = new QTableWidgetItem(text);
      item->setTextAlignment(Qt::AlignCenter);
      return item;
    };
    // construct the table
    ui_->table->setEditTriggers(QAbstractItemView::NoEditTriggers); // disable editing
    ui_->table->setRowCount(1);
    ui_->table->setColumnCount(5);
    ui_->table->setHorizontalHeaderLabels(QStringList() << "# Proteins" << "# Peptides" << "# Transitions" << "Decoy Frequency (%)" << "Reference Status");
    ui_->table->setItem(0, 0, getItem(QString::number(stats.protein_count)));
    ui_->table->setItem(0, 1, getItem(QString::number(stats.peptide_count)));
    ui_->table->setItem(0, 2, getItem(QString::number(stats.transition_count)));

    using TYPE = ReactionMonitoringTransition::DecoyTransitionType;
    auto count_copy = stats.decoy_counts; // allow to default construct missing values with 0 counts
    size_t all = count_copy[TYPE::DECOY] +
                 count_copy[TYPE::TARGET] +
                 count_copy[TYPE::UNKNOWN];
    if (all == 0)
    {
      all = 1; // avoid division by zero below
    }
    ui_->table->setItem(0, 3, getItem(QString::number(count_copy[TYPE::DECOY] * 100 / all)));
    ui_->table->setItem(0, 4, getItem((!stats.contains_invalid_references ? "valid" : "invalid")));
  }

  void SwathLibraryStats::updateFromFile(const QString& pqp_file)
  {
    TargetedExperiment te;
    TransitionPQPFile tr_file;
    tr_file.setLogType(ProgressLogger::GUI);
    tr_file.convertPQPToTargetedExperiment(pqp_file.toStdString().c_str(), te, true);
    //OpenSwath::LightTargetedExperiment transition_exp;
    //OpenSwathDataAccessHelper::convertTargetedExp(te, transition_exp);

    update(te.getSummary());
  }

} //namespace OpenMS



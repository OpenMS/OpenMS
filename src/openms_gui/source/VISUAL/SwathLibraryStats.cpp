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
    if (all == 0) all = 1; // avoid division by zero below
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



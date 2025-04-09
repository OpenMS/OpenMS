// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Dhanmoni Nath, Julianus Pfeuffer $
// --------------------------------------------------------------------------

#ifdef QT_WEBENGINEWIDGETS_LIB
#include <OpenMS/VISUAL/SequenceVisualizer.h>
#include <ui_SequenceVisualizer.h>

#include <QWebChannel>
#include <QString>

#include <QtWebEngineWidgets/QWebEngineView>

// This is the window that appears when we click on 'show' in the 'sequence' column of the protein table

namespace OpenMS
{
  SequenceVisualizer::SequenceVisualizer(QWidget* parent) :
      QWidget(parent), ui_(new Ui::SequenceVisualizer)
  {
    ui_->setupUi(this);
    view_ = new QWebEngineView(this);
    channel_ = new QWebChannel(&backend_); // setup Qt WebChannel API
    view_->page()->setWebChannel(channel_);
    channel_->registerObject(QString("Backend"), &backend_); // This object will be available in HTML file.
    view_->load(QUrl("qrc:/new/sequence_viz.html"));
    ui_->gridLayout->addWidget(view_);
  }

  SequenceVisualizer::~SequenceVisualizer()
  {
    channel_->deleteLater();
    view_->close();
    view_->deleteLater();
    delete ui_;
    deleteLater();
  }

  // Get protein and peptide data from the protein table and store inside the m_json_data_obj_ object. 
  // Inside the HTML file, this QObject will be available and we'll access these protein and 
  // peptide data using the qtWebEngine and webChannel API.
  void SequenceVisualizer::setProteinPeptideDataToJsonObj(
      const QString& accession_num,
      const QString& pro_seq,
      const QJsonArray& pep_data)
  {
    QJsonObject j;
    j["accession_num"] = accession_num;
    j["protein_sequence_data"] = pro_seq;
    j["peptides_data"] = pep_data;
    backend_.m_json_data_obj_ = std::move(j);
  }
}// namespace OpenMS
#endif
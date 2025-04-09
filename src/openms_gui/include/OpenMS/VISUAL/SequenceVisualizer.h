// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Dhanmoni Nath, Julianus Pfeuffer $
// --------------------------------------------------------------------------

#ifdef QT_WEBENGINEWIDGETS_LIB
#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>
#include <QJsonObject>

class QWebEngineView;
class QWebChannel;

namespace Ui
{
  class SequenceVisualizer;
}

namespace OpenMS
{
  class OPENMS_GUI_DLLAPI Backend : public QObject
  {
    Q_OBJECT

    // We can access the protein and peptide data using SequenceVisualizer.json_data_obj inside JS/HTML resource file
    Q_PROPERTY(QJsonObject json_data_obj MEMBER m_json_data_obj_ NOTIFY dataChanged_)
    signals:
      void dataChanged_();

  public:
    QJsonObject m_json_data_obj_;
  };

  class OPENMS_GUI_DLLAPI SequenceVisualizer : public QWidget
  {
    Q_OBJECT

  public:
    explicit SequenceVisualizer(QWidget* parent = nullptr);
    ~SequenceVisualizer() override;


  public slots:
    // this method sets protein and peptide data to m_json_data_obj_.
    void setProteinPeptideDataToJsonObj(
        const QString& accession_num, 
        const QString& pro_seq, 
        const QJsonArray& peptides_data);

  private:

    Ui::SequenceVisualizer* ui_;
    Backend backend_;
    QWebEngineView* view_;
    QWebChannel* channel_;
  };
}// namespace OpenMS
#endif
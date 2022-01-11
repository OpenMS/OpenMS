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
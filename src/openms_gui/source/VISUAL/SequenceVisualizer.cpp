#include <OpenMS/VISUAL/SequenceVisualizer.h>
#include <ui_SequenceVisualizer.h>

#include <QLabel>
#include <QDebug>
#include <QWebChannel>
#include <QPushButton>
#include <QMessageBox>
#include <QJsonArray>

#include <QtWebEngineWidgets/QWebEngineView>

namespace OpenMS
{
  SequenceVisualizer::SequenceVisualizer(QWidget* parent) :
      QWidget(parent), ui(new Ui::SequenceVisualizer)
  {
    ui->setupUi(this);
    QWebEngineView* view = new QWebEngineView(parent);
    QWebChannel* channel = new QWebChannel(this);
    view->page()->setWebChannel(channel);
    channel->registerObject(QString("SequenceVisualizer"), this);
    view->load(QUrl("qrc:/new/sequence_viz.html"));
    ui->gridLayout->addWidget(view);
  }
  SequenceVisualizer::~SequenceVisualizer()
  {
    delete ui;
  }
  void SequenceVisualizer::setProteinPeptideDataToJsonObj(QString pro_seq, QJsonArray pep_start_end_pos)
  {
    m_json_data_obj["protein_sequence_data"] = pro_seq;
    m_json_data_obj["peptides_start_end_pos"] = pep_start_end_pos;
   

    qDebug() << "sequence clicked";
  }

}// namespace OpenMS
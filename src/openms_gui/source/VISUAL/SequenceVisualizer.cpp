#include <OpenMS/VISUAL/SequenceVisualizer.h>
#include <ui_SequenceVisualizer.h>

#include <QLabel>
#include <QDebug>
#include <QWebChannel>
#include <QPushButton>
#include <QMessageBox>
#include <QtWebEngineWidgets/QWebEngineView>

namespace OpenMS
{
  SequenceVisualizer::SequenceVisualizer(QWidget* parent) :
      QWidget(parent), ui(new Ui::SequenceVisualizer)
  {
    ui->setupUi(this);
    QPushButton* btn = new QPushButton("Click here to visualize...");
    QWebEngineView* view = new QWebEngineView(parent);
    QWebChannel* channel = new QWebChannel(this);
    view->page()->setWebChannel(channel);
    channel->registerObject(QString("SequenceVisualizer"), this);
    connect(btn, SIGNAL(clicked()), this, SLOT(slotSendDataToJS()));
    view->load(QUrl("qrc:/new/sequence_viz.html"));
    
    ui->gridLayout->addWidget(btn);
    ui->gridLayout->addWidget(view);
  }
  SequenceVisualizer::~SequenceVisualizer()
  {
    delete ui;
  }
  void SequenceVisualizer::getData(const QString& pep_seq,const QJsonArray& accessionArr,const QJsonArray& sequenceArr)
  {
    m_json_data_obj["protein_accession_data"] = accessionArr;
    m_json_data_obj["protein_sequence_data"] = sequenceArr;
    m_json_data_obj["pep_seq"] = pep_seq;
    qDebug() << "getdata clicked";
  }

  void SequenceVisualizer::jscallme(const QString& datafromjs)
  {
    QMessageBox::information(NULL, "jscallme", "I'm called by js!");
    qDebug() << "getting data.." << datafromjs;
  }

  void SequenceVisualizer::slotSendDataToJS()
  {
    
   
  }
}// namespace OpenMS
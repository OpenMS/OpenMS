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
    /*bool res = connect(view, &QWebEngineView::loadFinished, this, &SequenceVisualizer::slotSendDataToJS);
    qDebug() << res << "!!<->!!";
    connect(view, &QObject::destroyed,
            [] { qDebug() << "Sender got deleted!"; });
    connect(this, &QObject::destroyed,
            [] { qDebug() << "Receiver got deleted!"; });*/
    ui->gridLayout->addWidget(btn);
    ui->gridLayout->addWidget(view);
  }
  SequenceVisualizer::~SequenceVisualizer()
  {
    delete ui;
  }
  void SequenceVisualizer::getData(QString seq)
  {
    m_sequence = seq;
    m_num = 20;
    qDebug() << "getdata clicked" << m_num << "<->" << m_sequence;
  }

  void SequenceVisualizer::jscallme(const QString& datafromjs)
  {
    QMessageBox::information(NULL, "jscallme", "I'm called by js!");
    qDebug() << "getting data.." << datafromjs;
  }

  void SequenceVisualizer::slotSendDataToJS()
  {
    
    qDebug() << "emitting signal..." << m_sequence;
    emit sendDataToJS(m_sequence);
  }
}// namespace OpenMS
#include <OpenMS/VISUAL/GoogleMap.h>
#include <ui_GoogleMap.h>

#include <QtWebEngineWidgets/QWebEngineView>
#include <QtWidgets/QPushButton>

namespace OpenMS
{
  GoogleMap::GoogleMap(QWidget* parent) :
      QWidget(parent), ui(new Ui::GoogleMap)
  {
    ui->setupUi(this);
    QWebEngineView* view = new QWebEngineView(parent);
    view->load(QUrl("qrc:/new/prefix/google_map"));
    QPushButton* button1 = new QPushButton("One");
    ui->verticalLayout->addWidget(view);
    ui->verticalLayout->addWidget(button1);
  }
  GoogleMap::~GoogleMap()
  {
    delete ui;
  }

}// namespace OpenMS
#include <OpenMS/VISUAL/PlotlyGraph.h>
#include <ui_PlotlyGraph.h>

#include <QtWebEngineWidgets/QWebEngineView>
#include <QtWidgets/QPushButton>

namespace OpenMS
{
  PlotlyGraph::PlotlyGraph(QWidget* parent) :
      QWidget(parent), ui(new Ui::PlotlyGraph)
  {
    ui->setupUi(this);
    QWebEngineView* view = new QWebEngineView(parent);
    view->load(QUrl("qrc:/new/prefix/plotly_graph"));
    ui->gridLayout->addWidget(view);
  }
  PlotlyGraph::~PlotlyGraph()
  {
    delete ui;
  }

}// namespace OpenMS
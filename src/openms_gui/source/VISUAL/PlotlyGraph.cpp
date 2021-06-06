#include <OpenMS/VISUAL/PlotlyGraph.h>
#include <ui_PlotlyGraph.h>

#include <QtWebEngineWidgets/QWebEngineView>
#include <QtWidgets/QPushButton>
#include <QWebChannel>
#include <QDebug>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

namespace OpenMS
{
  PlotlyGraph::PlotlyGraph(QWidget* parent) :
      QWidget(parent), ui(new Ui::PlotlyGraph)
  {
    ui->setupUi(this);
    QWebEngineView* view = new QWebEngineView(parent);
    QWebChannel* channel = new QWebChannel(this);
    view->page()->setWebChannel(channel);
    channel->registerObject(QString("plotlyGraph"), this);
    view->load(QUrl("qrc:/new/prefix/plotly_graph"));
    ui->gridLayout->addWidget(view);
  }
  PlotlyGraph::~PlotlyGraph()
  {
    delete ui;
  }

  void PlotlyGraph::on_AddLineBtn_clicked()
  {
    
    QJsonArray x_val;
    QJsonArray y_val;

    qDebug() << "add line button clicked!" << ui->xValues->text();
    QStringList xData = ui->xValues->text().split(",");
    QStringList yData = ui->yValues->text().split(",");

    foreach (QString num, xData)
    {
      x_val.push_back(num.toInt());
    }
    foreach (QString num, yData)
    {
      y_val.push_back(num.toInt());
    }

    if (x_val.size() == y_val.size())
    {
      QJsonObject xy_vals;
      xy_vals["x_val"] = x_val;
      xy_vals["y_val"] = y_val;

      QJsonObject obj;
      obj["plottingData"] = xy_vals;
      qDebug() << "correct number of input... Plotting graph :)";

      emit dataChanged(obj);
    }
    else
    {
      qDebug() << "incorrect number of input :(";
      QJsonObject err_obj;
      err_obj["error"] = "Please enter same number of x & y values";

      emit dataChanged(err_obj);
    }
  }

}// namespace OpenMS
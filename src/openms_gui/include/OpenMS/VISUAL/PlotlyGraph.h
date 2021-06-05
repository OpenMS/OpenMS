#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class PlotlyGraph;
}

namespace OpenMS
{

  class OPENMS_GUI_DLLAPI PlotlyGraph : public QWidget
  {
    Q_OBJECT
    //Q_PROPERTY(QString data MEMBER m_data NOTIFY dataChanged)
    //Q_PROPERTY(int num MEMBER m_num NOTIFY dataChanged)

  public:
    PlotlyGraph(QWidget* parent = nullptr);
    ~PlotlyGraph();
  public slots:
    //void jscallme(const QString& datafromjs);
    void on_AddLineBtn_clicked();
  signals:
    void dataChanged(const QJsonObject& val);


  private:
    Ui::PlotlyGraph* ui;
    //QString m_data;
    //int m_num;
  };
}// namespace OpenMS
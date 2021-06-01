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

  public:
    PlotlyGraph(QWidget* parent = nullptr);
    ~PlotlyGraph();


  private:
    Ui::PlotlyGraph* ui;
  };
}// namespace OpenMS
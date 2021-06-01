#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class GoogleMap;
}

namespace OpenMS
{

  class OPENMS_GUI_DLLAPI GoogleMap : public QWidget
  {
    Q_OBJECT

  public:
    GoogleMap(QWidget* parent = nullptr);
    ~GoogleMap();


  private:
    Ui::GoogleMap* ui;
  };
}// namespace OpenMS
#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class FancyIcon;
}

namespace OpenMS
{

  class OPENMS_GUI_DLLAPI FancyIcon : public QWidget
  {
    Q_OBJECT

  public:
    FancyIcon(QWidget* parent = nullptr);
    ~FancyIcon();


  private:
    Ui::FancyIcon* ui;
  };
}// namespace OpenMS
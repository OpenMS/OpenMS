#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class Editor;
}

namespace OpenMS
{
  class OPENMS_GUI_DLLAPI Editor : public QWidget
  {
    Q_OBJECT

  public:
    Editor(QWidget* parent = nullptr);
    ~Editor();

  private slots:
    void on_New_clicked();

    void on_Open_clicked();

    void on_Save_clicked();

    void on_Save_As_clicked();

    void on_Copy_clicked();

    void on_Paste_clicked();

    void on_Cut_clicked();

    void on_Undo_clicked();

    void on_Redo_clicked();
  private:
    Ui::Editor* ui;
    QString file_path_;
  };
}
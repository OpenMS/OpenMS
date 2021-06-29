#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class SequenceVisualizer;
}

namespace OpenMS
{

  class OPENMS_GUI_DLLAPI SequenceVisualizer : public QWidget
  {
    Q_OBJECT
    Q_PROPERTY(QString sequence MEMBER m_sequence)
    Q_PROPERTY(int num MEMBER m_num)

  public:
    SequenceVisualizer(QWidget* parent = nullptr);
    ~SequenceVisualizer();

  public slots:

    void jscallme(const QString& datafromjs);
    void slotSendDataToJS();
    void getData(QString txt);

  signals:
    void sendDataToJS(const QString&);

  private:
    Ui::SequenceVisualizer* ui;
    QString m_sequence;
    int m_num = 10;
  };
}// namespace OpenMS
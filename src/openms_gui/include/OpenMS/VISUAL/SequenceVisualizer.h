#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>
#include <QJsonObject>

namespace Ui
{
  class SequenceVisualizer;
}

namespace OpenMS
{

  class OPENMS_GUI_DLLAPI SequenceVisualizer : public QWidget
  {
    Q_OBJECT
    Q_PROPERTY(QJsonObject json_data_obj MEMBER m_json_data_obj)

  public:
    SequenceVisualizer(QWidget* parent = nullptr);
    ~SequenceVisualizer();

  public slots:

    void jscallme(const QString& datafromjs);
    void slotSendDataToJS();
    void getData(const QString& txt,const QJsonArray& accessionArr,const QJsonArray& sequenceArr);

  signals:
    void sendDataToJS(const QString&);

  private:
    Ui::SequenceVisualizer* ui;
    QJsonObject m_json_data_obj;
  };
}// namespace OpenMS
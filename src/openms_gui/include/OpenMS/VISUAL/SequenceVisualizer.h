#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>
#include <QJsonObject>

#include <unordered_map>
#include <vector>

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

    void setProteinPeptideDataToJsonObj(QString pro_seq, QJsonArray peptides_data, QJsonArray pep_mod_data);

  private:
    Ui::SequenceVisualizer* ui;
    QJsonObject m_json_data_obj;
  };
}// namespace OpenMS
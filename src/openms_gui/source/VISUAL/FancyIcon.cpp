#include <OpenMS/VISUAL/FancyIcon.h>
#include <ui_FancyIcon.h>
#include <QIcon>

namespace OpenMS 
{
  FancyIcon::FancyIcon(QWidget* parent) :
      QWidget(parent), ui(new Ui::FancyIcon)
  {
    ui->setupUi(this);
    ui->pushButton->setIcon(QIcon(":/new/TOPPView.png"));
  }
  FancyIcon::~FancyIcon()
  {
    delete ui;
  }

}
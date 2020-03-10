#include "SwathTabWidget.h"
#include "ui_SwathTabWidget.h"

SwathTabWidget::SwathTabWidget(QWidget *parent) :
    QTabWidget(parent),
    ui(new Ui::SwathTabWidget)
{
    ui->setupUi(this);
}

SwathTabWidget::~SwathTabWidget()
{
    delete ui;
}

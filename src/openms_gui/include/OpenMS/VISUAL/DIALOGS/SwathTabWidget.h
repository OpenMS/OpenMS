#ifndef SWATHTABWIDGET_H
#define SWATHTABWIDGET_H

#include <QTabWidget>

namespace Ui {
class SwathTabWidget;
}

class SwathTabWidget : public QTabWidget
{
    Q_OBJECT

public:
    explicit SwathTabWidget(QWidget *parent = nullptr);
    ~SwathTabWidget();

private:
    Ui::SwathTabWidget *ui;
};

#endif // SWATHTABWIDGET_H

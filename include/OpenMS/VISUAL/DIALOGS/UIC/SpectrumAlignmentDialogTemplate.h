/********************************************************************************
** Form generated from reading ui file 'SpectrumAlignmentDialog.ui'
**
** Created: Tue Dec 9 19:34:18 2008
**      by: Qt User Interface Compiler version 4.3.5
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef SPECTRUMALIGNMENTDIALOGTEMPLATE_H
#define SPECTRUMALIGNMENTDIALOGTEMPLATE_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QLabel>
#include <QtGui/QRadioButton>
#include <QtGui/QVBoxLayout>

class Ui_SpectrumAlignmentDialogTemplate
{
public:
    QVBoxLayout *vboxLayout;
    QHBoxLayout *hboxLayout;
    QLabel *tolerance_label;
    QDoubleSpinBox *tolerance_spinbox;
    QGroupBox *unit_group;
    QVBoxLayout *vboxLayout1;
    QRadioButton *da;
    QRadioButton *ppm;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *SpectrumAlignmentDialogTemplate)
    {
    if (SpectrumAlignmentDialogTemplate->objectName().isEmpty())
        SpectrumAlignmentDialogTemplate->setObjectName(QString::fromUtf8("SpectrumAlignmentDialogTemplate"));
    SpectrumAlignmentDialogTemplate->resize(178, 170);
    vboxLayout = new QVBoxLayout(SpectrumAlignmentDialogTemplate);
    vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
    hboxLayout = new QHBoxLayout();
    hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
    tolerance_label = new QLabel(SpectrumAlignmentDialogTemplate);
    tolerance_label->setObjectName(QString::fromUtf8("tolerance_label"));

    hboxLayout->addWidget(tolerance_label);

    tolerance_spinbox = new QDoubleSpinBox(SpectrumAlignmentDialogTemplate);
    tolerance_spinbox->setObjectName(QString::fromUtf8("tolerance_spinbox"));
    tolerance_spinbox->setMaximum(1e+06);
    tolerance_spinbox->setSingleStep(0.01);
    tolerance_spinbox->setValue(0.3);

    hboxLayout->addWidget(tolerance_spinbox);


    vboxLayout->addLayout(hboxLayout);

    unit_group = new QGroupBox(SpectrumAlignmentDialogTemplate);
    unit_group->setObjectName(QString::fromUtf8("unit_group"));
    vboxLayout1 = new QVBoxLayout(unit_group);
    vboxLayout1->setObjectName(QString::fromUtf8("vboxLayout1"));
    da = new QRadioButton(unit_group);
    da->setObjectName(QString::fromUtf8("da"));

    vboxLayout1->addWidget(da);

    ppm = new QRadioButton(unit_group);
    ppm->setObjectName(QString::fromUtf8("ppm"));

    vboxLayout1->addWidget(ppm);


    vboxLayout->addWidget(unit_group);

    buttonBox = new QDialogButtonBox(SpectrumAlignmentDialogTemplate);
    buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
    buttonBox->setOrientation(Qt::Horizontal);
    buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

    vboxLayout->addWidget(buttonBox);


    retranslateUi(SpectrumAlignmentDialogTemplate);
    QObject::connect(buttonBox, SIGNAL(accepted()), SpectrumAlignmentDialogTemplate, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), SpectrumAlignmentDialogTemplate, SLOT(reject()));

    QMetaObject::connectSlotsByName(SpectrumAlignmentDialogTemplate);
    } // setupUi

    void retranslateUi(QDialog *SpectrumAlignmentDialogTemplate)
    {
    SpectrumAlignmentDialogTemplate->setWindowTitle(QApplication::translate("SpectrumAlignmentDialogTemplate", "Spectrum alignment", 0, QApplication::UnicodeUTF8));
    tolerance_label->setText(QApplication::translate("SpectrumAlignmentDialogTemplate", "Tolerance:", 0, QApplication::UnicodeUTF8));
    unit_group->setTitle(QApplication::translate("SpectrumAlignmentDialogTemplate", "Unit", 0, QApplication::UnicodeUTF8));
    da->setText(QApplication::translate("SpectrumAlignmentDialogTemplate", "Da", 0, QApplication::UnicodeUTF8));
    ppm->setText(QApplication::translate("SpectrumAlignmentDialogTemplate", "ppm", 0, QApplication::UnicodeUTF8));
    Q_UNUSED(SpectrumAlignmentDialogTemplate);
    } // retranslateUi

};

namespace Ui {
    class SpectrumAlignmentDialogTemplate: public Ui_SpectrumAlignmentDialogTemplate {};
} // namespace Ui

#endif // SPECTRUMALIGNMENTDIALOGTEMPLATE_H

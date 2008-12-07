/********************************************************************************
** Form generated from reading ui file 'SpectrumAlignmentDialog.ui'
**
** Created: Sun Dec 7 22:10:29 2008
**      by: Qt User Interface Compiler version 4.4.3
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

QT_BEGIN_NAMESPACE

class Ui_SpectrumAlignmentDialogTemplate
{
public:
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout;
    QLabel *tolerance_label;
    QDoubleSpinBox *tolerance_spinbox;
    QGroupBox *unit_group;
    QVBoxLayout *verticalLayout;
    QRadioButton *da;
    QRadioButton *ppm;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *SpectrumAlignmentDialogTemplate)
    {
    if (SpectrumAlignmentDialogTemplate->objectName().isEmpty())
        SpectrumAlignmentDialogTemplate->setObjectName(QString::fromUtf8("SpectrumAlignmentDialogTemplate"));
    SpectrumAlignmentDialogTemplate->resize(178, 170);
    verticalLayout_2 = new QVBoxLayout(SpectrumAlignmentDialogTemplate);
    verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
    horizontalLayout = new QHBoxLayout();
    horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
    tolerance_label = new QLabel(SpectrumAlignmentDialogTemplate);
    tolerance_label->setObjectName(QString::fromUtf8("tolerance_label"));

    horizontalLayout->addWidget(tolerance_label);

    tolerance_spinbox = new QDoubleSpinBox(SpectrumAlignmentDialogTemplate);
    tolerance_spinbox->setObjectName(QString::fromUtf8("tolerance_spinbox"));
    tolerance_spinbox->setMaximum(1e+06);
    tolerance_spinbox->setSingleStep(0.01);
    tolerance_spinbox->setValue(0.3);

    horizontalLayout->addWidget(tolerance_spinbox);


    verticalLayout_2->addLayout(horizontalLayout);

    unit_group = new QGroupBox(SpectrumAlignmentDialogTemplate);
    unit_group->setObjectName(QString::fromUtf8("unit_group"));
    verticalLayout = new QVBoxLayout(unit_group);
    verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
    da = new QRadioButton(unit_group);
    da->setObjectName(QString::fromUtf8("da"));

    verticalLayout->addWidget(da);

    ppm = new QRadioButton(unit_group);
    ppm->setObjectName(QString::fromUtf8("ppm"));

    verticalLayout->addWidget(ppm);


    verticalLayout_2->addWidget(unit_group);

    buttonBox = new QDialogButtonBox(SpectrumAlignmentDialogTemplate);
    buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
    buttonBox->setOrientation(Qt::Horizontal);
    buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

    verticalLayout_2->addWidget(buttonBox);


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

QT_END_NAMESPACE

#endif // SPECTRUMALIGNMENTDIALOGTEMPLATE_H

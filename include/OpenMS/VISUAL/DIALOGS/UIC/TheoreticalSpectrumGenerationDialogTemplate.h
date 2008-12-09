/********************************************************************************
** Form generated from reading ui file 'TheoreticalSpectrumGenerationDialog.ui'
**
** Created: Tue Dec 9 19:34:18 2008
**      by: Qt User Interface Compiler version 4.3.5
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef THEORETICALSPECTRUMGENERATIONDIALOGTEMPLATE_H
#define THEORETICALSPECTRUMGENERATIONDIALOGTEMPLATE_H

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
#include <QtGui/QLineEdit>
#include <QtGui/QListWidget>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>

class Ui_TheoreticalSpectrumGenerationDialogTemplate
{
public:
    QVBoxLayout *vboxLayout;
    QVBoxLayout *vboxLayout1;
    QLabel *enter_seq_label;
    QLineEdit *line_edit;
    QHBoxLayout *hboxLayout;
    QVBoxLayout *vboxLayout2;
    QHBoxLayout *hboxLayout1;
    QHBoxLayout *hboxLayout2;
    QLabel *spin_label;
    QSpinBox *spin_box;
    QHBoxLayout *hboxLayout3;
    QLabel *max_iso_label;
    QSpinBox *max_iso_spinbox;
    QVBoxLayout *vboxLayout3;
    QLabel *generate_label;
    QListWidget *list_widget;
    QGroupBox *groupBox;
    QVBoxLayout *vboxLayout4;
    QHBoxLayout *hboxLayout4;
    QLabel *label;
    QDoubleSpinBox *a_intensity;
    QHBoxLayout *hboxLayout5;
    QLabel *label_2;
    QDoubleSpinBox *b_intensity;
    QHBoxLayout *hboxLayout6;
    QLabel *label_3;
    QDoubleSpinBox *c_intensity;
    QHBoxLayout *hboxLayout7;
    QLabel *label_4;
    QDoubleSpinBox *x_intensity;
    QHBoxLayout *hboxLayout8;
    QLabel *label_5;
    QDoubleSpinBox *y_intensity;
    QHBoxLayout *hboxLayout9;
    QLabel *label_6;
    QDoubleSpinBox *z_intensity;
    QHBoxLayout *hboxLayout10;
    QLabel *label_7;
    QSpinBox *rel_loss_intensity;
    QDialogButtonBox *button_box;

    void setupUi(QDialog *TheoreticalSpectrumGenerationDialogTemplate)
    {
    if (TheoreticalSpectrumGenerationDialogTemplate->objectName().isEmpty())
        TheoreticalSpectrumGenerationDialogTemplate->setObjectName(QString::fromUtf8("TheoreticalSpectrumGenerationDialogTemplate"));
    TheoreticalSpectrumGenerationDialogTemplate->resize(450, 341);
    vboxLayout = new QVBoxLayout(TheoreticalSpectrumGenerationDialogTemplate);
    vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
    vboxLayout1 = new QVBoxLayout();
    vboxLayout1->setObjectName(QString::fromUtf8("vboxLayout1"));
    enter_seq_label = new QLabel(TheoreticalSpectrumGenerationDialogTemplate);
    enter_seq_label->setObjectName(QString::fromUtf8("enter_seq_label"));

    vboxLayout1->addWidget(enter_seq_label);

    line_edit = new QLineEdit(TheoreticalSpectrumGenerationDialogTemplate);
    line_edit->setObjectName(QString::fromUtf8("line_edit"));
    QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(line_edit->sizePolicy().hasHeightForWidth());
    line_edit->setSizePolicy(sizePolicy);

    vboxLayout1->addWidget(line_edit);


    vboxLayout->addLayout(vboxLayout1);

    hboxLayout = new QHBoxLayout();
    hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
    vboxLayout2 = new QVBoxLayout();
    vboxLayout2->setObjectName(QString::fromUtf8("vboxLayout2"));
    hboxLayout1 = new QHBoxLayout();
    hboxLayout1->setObjectName(QString::fromUtf8("hboxLayout1"));
    hboxLayout2 = new QHBoxLayout();
    hboxLayout2->setObjectName(QString::fromUtf8("hboxLayout2"));
    spin_label = new QLabel(TheoreticalSpectrumGenerationDialogTemplate);
    spin_label->setObjectName(QString::fromUtf8("spin_label"));
    QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Fixed);
    sizePolicy1.setHorizontalStretch(0);
    sizePolicy1.setVerticalStretch(0);
    sizePolicy1.setHeightForWidth(spin_label->sizePolicy().hasHeightForWidth());
    spin_label->setSizePolicy(sizePolicy1);

    hboxLayout2->addWidget(spin_label);

    spin_box = new QSpinBox(TheoreticalSpectrumGenerationDialogTemplate);
    spin_box->setObjectName(QString::fromUtf8("spin_box"));
    QSizePolicy sizePolicy2(QSizePolicy::Maximum, QSizePolicy::Fixed);
    sizePolicy2.setHorizontalStretch(0);
    sizePolicy2.setVerticalStretch(0);
    sizePolicy2.setHeightForWidth(spin_box->sizePolicy().hasHeightForWidth());
    spin_box->setSizePolicy(sizePolicy2);
    spin_box->setValue(1);

    hboxLayout2->addWidget(spin_box);


    hboxLayout1->addLayout(hboxLayout2);

    hboxLayout3 = new QHBoxLayout();
    hboxLayout3->setObjectName(QString::fromUtf8("hboxLayout3"));
    max_iso_label = new QLabel(TheoreticalSpectrumGenerationDialogTemplate);
    max_iso_label->setObjectName(QString::fromUtf8("max_iso_label"));
    max_iso_label->setEnabled(false);
    max_iso_label->setMaximumSize(QSize(16777213, 16777215));

    hboxLayout3->addWidget(max_iso_label);

    max_iso_spinbox = new QSpinBox(TheoreticalSpectrumGenerationDialogTemplate);
    max_iso_spinbox->setObjectName(QString::fromUtf8("max_iso_spinbox"));
    max_iso_spinbox->setEnabled(false);
    max_iso_spinbox->setValue(2);

    hboxLayout3->addWidget(max_iso_spinbox);


    hboxLayout1->addLayout(hboxLayout3);


    vboxLayout2->addLayout(hboxLayout1);

    vboxLayout3 = new QVBoxLayout();
    vboxLayout3->setObjectName(QString::fromUtf8("vboxLayout3"));
    generate_label = new QLabel(TheoreticalSpectrumGenerationDialogTemplate);
    generate_label->setObjectName(QString::fromUtf8("generate_label"));
    sizePolicy1.setHeightForWidth(generate_label->sizePolicy().hasHeightForWidth());
    generate_label->setSizePolicy(sizePolicy1);

    vboxLayout3->addWidget(generate_label);

    list_widget = new QListWidget(TheoreticalSpectrumGenerationDialogTemplate);
    list_widget->setObjectName(QString::fromUtf8("list_widget"));
    QSizePolicy sizePolicy3(QSizePolicy::Ignored, QSizePolicy::Ignored);
    sizePolicy3.setHorizontalStretch(0);
    sizePolicy3.setVerticalStretch(0);
    sizePolicy3.setHeightForWidth(list_widget->sizePolicy().hasHeightForWidth());
    list_widget->setSizePolicy(sizePolicy3);
    list_widget->setMinimumSize(QSize(100, 155));
    list_widget->setSelectionMode(QAbstractItemView::NoSelection);

    vboxLayout3->addWidget(list_widget);


    vboxLayout2->addLayout(vboxLayout3);


    hboxLayout->addLayout(vboxLayout2);

    groupBox = new QGroupBox(TheoreticalSpectrumGenerationDialogTemplate);
    groupBox->setObjectName(QString::fromUtf8("groupBox"));
    vboxLayout4 = new QVBoxLayout(groupBox);
    vboxLayout4->setObjectName(QString::fromUtf8("vboxLayout4"));
    hboxLayout4 = new QHBoxLayout();
    hboxLayout4->setObjectName(QString::fromUtf8("hboxLayout4"));
    label = new QLabel(groupBox);
    label->setObjectName(QString::fromUtf8("label"));

    hboxLayout4->addWidget(label);

    a_intensity = new QDoubleSpinBox(groupBox);
    a_intensity->setObjectName(QString::fromUtf8("a_intensity"));
    a_intensity->setValue(1);

    hboxLayout4->addWidget(a_intensity);


    vboxLayout4->addLayout(hboxLayout4);

    hboxLayout5 = new QHBoxLayout();
    hboxLayout5->setObjectName(QString::fromUtf8("hboxLayout5"));
    label_2 = new QLabel(groupBox);
    label_2->setObjectName(QString::fromUtf8("label_2"));

    hboxLayout5->addWidget(label_2);

    b_intensity = new QDoubleSpinBox(groupBox);
    b_intensity->setObjectName(QString::fromUtf8("b_intensity"));
    b_intensity->setValue(1);

    hboxLayout5->addWidget(b_intensity);


    vboxLayout4->addLayout(hboxLayout5);

    hboxLayout6 = new QHBoxLayout();
    hboxLayout6->setObjectName(QString::fromUtf8("hboxLayout6"));
    label_3 = new QLabel(groupBox);
    label_3->setObjectName(QString::fromUtf8("label_3"));

    hboxLayout6->addWidget(label_3);

    c_intensity = new QDoubleSpinBox(groupBox);
    c_intensity->setObjectName(QString::fromUtf8("c_intensity"));
    c_intensity->setValue(1);

    hboxLayout6->addWidget(c_intensity);


    vboxLayout4->addLayout(hboxLayout6);

    hboxLayout7 = new QHBoxLayout();
    hboxLayout7->setObjectName(QString::fromUtf8("hboxLayout7"));
    label_4 = new QLabel(groupBox);
    label_4->setObjectName(QString::fromUtf8("label_4"));

    hboxLayout7->addWidget(label_4);

    x_intensity = new QDoubleSpinBox(groupBox);
    x_intensity->setObjectName(QString::fromUtf8("x_intensity"));
    x_intensity->setValue(1);

    hboxLayout7->addWidget(x_intensity);


    vboxLayout4->addLayout(hboxLayout7);

    hboxLayout8 = new QHBoxLayout();
    hboxLayout8->setObjectName(QString::fromUtf8("hboxLayout8"));
    label_5 = new QLabel(groupBox);
    label_5->setObjectName(QString::fromUtf8("label_5"));

    hboxLayout8->addWidget(label_5);

    y_intensity = new QDoubleSpinBox(groupBox);
    y_intensity->setObjectName(QString::fromUtf8("y_intensity"));
    y_intensity->setValue(1);

    hboxLayout8->addWidget(y_intensity);


    vboxLayout4->addLayout(hboxLayout8);

    hboxLayout9 = new QHBoxLayout();
    hboxLayout9->setObjectName(QString::fromUtf8("hboxLayout9"));
    label_6 = new QLabel(groupBox);
    label_6->setObjectName(QString::fromUtf8("label_6"));

    hboxLayout9->addWidget(label_6);

    z_intensity = new QDoubleSpinBox(groupBox);
    z_intensity->setObjectName(QString::fromUtf8("z_intensity"));
    z_intensity->setValue(1);

    hboxLayout9->addWidget(z_intensity);


    vboxLayout4->addLayout(hboxLayout9);

    hboxLayout10 = new QHBoxLayout();
    hboxLayout10->setObjectName(QString::fromUtf8("hboxLayout10"));
    label_7 = new QLabel(groupBox);
    label_7->setObjectName(QString::fromUtf8("label_7"));

    hboxLayout10->addWidget(label_7);

    rel_loss_intensity = new QSpinBox(groupBox);
    rel_loss_intensity->setObjectName(QString::fromUtf8("rel_loss_intensity"));
    rel_loss_intensity->setMaximum(100);
    rel_loss_intensity->setValue(10);

    hboxLayout10->addWidget(rel_loss_intensity);


    vboxLayout4->addLayout(hboxLayout10);


    hboxLayout->addWidget(groupBox);


    vboxLayout->addLayout(hboxLayout);

    button_box = new QDialogButtonBox(TheoreticalSpectrumGenerationDialogTemplate);
    button_box->setObjectName(QString::fromUtf8("button_box"));
    QSizePolicy sizePolicy4(QSizePolicy::Ignored, QSizePolicy::Fixed);
    sizePolicy4.setHorizontalStretch(0);
    sizePolicy4.setVerticalStretch(0);
    sizePolicy4.setHeightForWidth(button_box->sizePolicy().hasHeightForWidth());
    button_box->setSizePolicy(sizePolicy4);
    button_box->setOrientation(Qt::Horizontal);
    button_box->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

    vboxLayout->addWidget(button_box);


    retranslateUi(TheoreticalSpectrumGenerationDialogTemplate);
    QObject::connect(button_box, SIGNAL(accepted()), TheoreticalSpectrumGenerationDialogTemplate, SLOT(accept()));
    QObject::connect(button_box, SIGNAL(rejected()), TheoreticalSpectrumGenerationDialogTemplate, SLOT(reject()));

    QMetaObject::connectSlotsByName(TheoreticalSpectrumGenerationDialogTemplate);
    } // setupUi

    void retranslateUi(QDialog *TheoreticalSpectrumGenerationDialogTemplate)
    {
    TheoreticalSpectrumGenerationDialogTemplate->setWindowTitle(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Generate theoretical spectrum", 0, QApplication::UnicodeUTF8));
    enter_seq_label->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Enter peptide sequence:", 0, QApplication::UnicodeUTF8));
    spin_label->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Charge:", 0, QApplication::UnicodeUTF8));
    max_iso_label->setWhatsThis(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Defines the maximal isotopic peak which is added", 0, QApplication::UnicodeUTF8));
    max_iso_label->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Max Isotope:", 0, QApplication::UnicodeUTF8));
    generate_label->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Generate:", 0, QApplication::UnicodeUTF8));
    list_widget->clear();

    QListWidgetItem *__item = new QListWidgetItem(list_widget);
    __item->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "A-ions", 0, QApplication::UnicodeUTF8));

    QListWidgetItem *__item1 = new QListWidgetItem(list_widget);
    __item1->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "B-ions", 0, QApplication::UnicodeUTF8));

    QListWidgetItem *__item2 = new QListWidgetItem(list_widget);
    __item2->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "C-ions", 0, QApplication::UnicodeUTF8));

    QListWidgetItem *__item3 = new QListWidgetItem(list_widget);
    __item3->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "X-ions", 0, QApplication::UnicodeUTF8));

    QListWidgetItem *__item4 = new QListWidgetItem(list_widget);
    __item4->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Y-ions", 0, QApplication::UnicodeUTF8));

    QListWidgetItem *__item5 = new QListWidgetItem(list_widget);
    __item5->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Z-ions", 0, QApplication::UnicodeUTF8));

    QListWidgetItem *__item6 = new QListWidgetItem(list_widget);
    __item6->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Precursor", 0, QApplication::UnicodeUTF8));

    QListWidgetItem *__item7 = new QListWidgetItem(list_widget);
    __item7->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Neutral losses", 0, QApplication::UnicodeUTF8));

    QListWidgetItem *__item8 = new QListWidgetItem(list_widget);
    __item8->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Isotope clusters", 0, QApplication::UnicodeUTF8));
    groupBox->setTitle(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Intensities", 0, QApplication::UnicodeUTF8));
    label->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "A-ions", 0, QApplication::UnicodeUTF8));
    label_2->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "B-ions", 0, QApplication::UnicodeUTF8));
    label_3->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "C-ions", 0, QApplication::UnicodeUTF8));
    label_4->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "X-ions", 0, QApplication::UnicodeUTF8));
    label_5->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Y-ions", 0, QApplication::UnicodeUTF8));
    label_6->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Z-ions", 0, QApplication::UnicodeUTF8));
    label_7->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Relative loss in %", 0, QApplication::UnicodeUTF8));
    Q_UNUSED(TheoreticalSpectrumGenerationDialogTemplate);
    } // retranslateUi

};

namespace Ui {
    class TheoreticalSpectrumGenerationDialogTemplate: public Ui_TheoreticalSpectrumGenerationDialogTemplate {};
} // namespace Ui

#endif // THEORETICALSPECTRUMGENERATIONDIALOGTEMPLATE_H

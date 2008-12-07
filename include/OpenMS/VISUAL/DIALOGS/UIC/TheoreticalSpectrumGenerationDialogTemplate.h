/********************************************************************************
** Form generated from reading ui file 'TheoreticalSpectrumGenerationDialog.ui'
**
** Created: Sun Dec 7 22:10:29 2008
**      by: Qt User Interface Compiler version 4.4.3
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

QT_BEGIN_NAMESPACE

class Ui_TheoreticalSpectrumGenerationDialogTemplate
{
public:
    QVBoxLayout *verticalLayout_5;
    QVBoxLayout *verticalLayout;
    QLabel *enter_seq_label;
    QLineEdit *line_edit;
    QHBoxLayout *horizontalLayout_11;
    QVBoxLayout *verticalLayout_4;
    QHBoxLayout *horizontalLayout_10;
    QHBoxLayout *horizontalLayout;
    QLabel *spin_label;
    QSpinBox *spin_box;
    QHBoxLayout *horizontalLayout_2;
    QLabel *max_iso_label;
    QSpinBox *max_iso_spinbox;
    QVBoxLayout *verticalLayout_2;
    QLabel *generate_label;
    QListWidget *list_widget;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout_3;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label;
    QDoubleSpinBox *a_intensity;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_2;
    QDoubleSpinBox *b_intensity;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_3;
    QDoubleSpinBox *c_intensity;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_4;
    QDoubleSpinBox *x_intensity;
    QHBoxLayout *horizontalLayout_7;
    QLabel *label_5;
    QDoubleSpinBox *y_intensity;
    QHBoxLayout *horizontalLayout_8;
    QLabel *label_6;
    QDoubleSpinBox *z_intensity;
    QHBoxLayout *horizontalLayout_9;
    QLabel *label_7;
    QSpinBox *rel_loss_intensity;
    QDialogButtonBox *button_box;

    void setupUi(QDialog *TheoreticalSpectrumGenerationDialogTemplate)
    {
    if (TheoreticalSpectrumGenerationDialogTemplate->objectName().isEmpty())
        TheoreticalSpectrumGenerationDialogTemplate->setObjectName(QString::fromUtf8("TheoreticalSpectrumGenerationDialogTemplate"));
    TheoreticalSpectrumGenerationDialogTemplate->resize(450, 341);
    verticalLayout_5 = new QVBoxLayout(TheoreticalSpectrumGenerationDialogTemplate);
    verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
    verticalLayout = new QVBoxLayout();
    verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
    enter_seq_label = new QLabel(TheoreticalSpectrumGenerationDialogTemplate);
    enter_seq_label->setObjectName(QString::fromUtf8("enter_seq_label"));

    verticalLayout->addWidget(enter_seq_label);

    line_edit = new QLineEdit(TheoreticalSpectrumGenerationDialogTemplate);
    line_edit->setObjectName(QString::fromUtf8("line_edit"));
    QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(line_edit->sizePolicy().hasHeightForWidth());
    line_edit->setSizePolicy(sizePolicy);

    verticalLayout->addWidget(line_edit);


    verticalLayout_5->addLayout(verticalLayout);

    horizontalLayout_11 = new QHBoxLayout();
    horizontalLayout_11->setObjectName(QString::fromUtf8("horizontalLayout_11"));
    verticalLayout_4 = new QVBoxLayout();
    verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
    horizontalLayout_10 = new QHBoxLayout();
    horizontalLayout_10->setObjectName(QString::fromUtf8("horizontalLayout_10"));
    horizontalLayout = new QHBoxLayout();
    horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
    spin_label = new QLabel(TheoreticalSpectrumGenerationDialogTemplate);
    spin_label->setObjectName(QString::fromUtf8("spin_label"));
    QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Fixed);
    sizePolicy1.setHorizontalStretch(0);
    sizePolicy1.setVerticalStretch(0);
    sizePolicy1.setHeightForWidth(spin_label->sizePolicy().hasHeightForWidth());
    spin_label->setSizePolicy(sizePolicy1);

    horizontalLayout->addWidget(spin_label);

    spin_box = new QSpinBox(TheoreticalSpectrumGenerationDialogTemplate);
    spin_box->setObjectName(QString::fromUtf8("spin_box"));
    QSizePolicy sizePolicy2(QSizePolicy::Maximum, QSizePolicy::Fixed);
    sizePolicy2.setHorizontalStretch(0);
    sizePolicy2.setVerticalStretch(0);
    sizePolicy2.setHeightForWidth(spin_box->sizePolicy().hasHeightForWidth());
    spin_box->setSizePolicy(sizePolicy2);
    spin_box->setValue(1);

    horizontalLayout->addWidget(spin_box);


    horizontalLayout_10->addLayout(horizontalLayout);

    horizontalLayout_2 = new QHBoxLayout();
    horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
    max_iso_label = new QLabel(TheoreticalSpectrumGenerationDialogTemplate);
    max_iso_label->setObjectName(QString::fromUtf8("max_iso_label"));
    max_iso_label->setEnabled(false);
    max_iso_label->setMaximumSize(QSize(16777213, 16777215));

    horizontalLayout_2->addWidget(max_iso_label);

    max_iso_spinbox = new QSpinBox(TheoreticalSpectrumGenerationDialogTemplate);
    max_iso_spinbox->setObjectName(QString::fromUtf8("max_iso_spinbox"));
    max_iso_spinbox->setEnabled(false);
    max_iso_spinbox->setValue(2);

    horizontalLayout_2->addWidget(max_iso_spinbox);


    horizontalLayout_10->addLayout(horizontalLayout_2);


    verticalLayout_4->addLayout(horizontalLayout_10);

    verticalLayout_2 = new QVBoxLayout();
    verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
    generate_label = new QLabel(TheoreticalSpectrumGenerationDialogTemplate);
    generate_label->setObjectName(QString::fromUtf8("generate_label"));
    sizePolicy1.setHeightForWidth(generate_label->sizePolicy().hasHeightForWidth());
    generate_label->setSizePolicy(sizePolicy1);

    verticalLayout_2->addWidget(generate_label);

    list_widget = new QListWidget(TheoreticalSpectrumGenerationDialogTemplate);
    new QListWidgetItem(list_widget);
    new QListWidgetItem(list_widget);
    new QListWidgetItem(list_widget);
    new QListWidgetItem(list_widget);
    new QListWidgetItem(list_widget);
    new QListWidgetItem(list_widget);
    new QListWidgetItem(list_widget);
    new QListWidgetItem(list_widget);
    new QListWidgetItem(list_widget);
    list_widget->setObjectName(QString::fromUtf8("list_widget"));
    QSizePolicy sizePolicy3(QSizePolicy::Ignored, QSizePolicy::Ignored);
    sizePolicy3.setHorizontalStretch(0);
    sizePolicy3.setVerticalStretch(0);
    sizePolicy3.setHeightForWidth(list_widget->sizePolicy().hasHeightForWidth());
    list_widget->setSizePolicy(sizePolicy3);
    list_widget->setMinimumSize(QSize(100, 155));
    list_widget->setSelectionMode(QAbstractItemView::NoSelection);

    verticalLayout_2->addWidget(list_widget);


    verticalLayout_4->addLayout(verticalLayout_2);


    horizontalLayout_11->addLayout(verticalLayout_4);

    groupBox = new QGroupBox(TheoreticalSpectrumGenerationDialogTemplate);
    groupBox->setObjectName(QString::fromUtf8("groupBox"));
    verticalLayout_3 = new QVBoxLayout(groupBox);
    verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
    horizontalLayout_3 = new QHBoxLayout();
    horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
    label = new QLabel(groupBox);
    label->setObjectName(QString::fromUtf8("label"));

    horizontalLayout_3->addWidget(label);

    a_intensity = new QDoubleSpinBox(groupBox);
    a_intensity->setObjectName(QString::fromUtf8("a_intensity"));
    a_intensity->setValue(1);

    horizontalLayout_3->addWidget(a_intensity);


    verticalLayout_3->addLayout(horizontalLayout_3);

    horizontalLayout_4 = new QHBoxLayout();
    horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
    label_2 = new QLabel(groupBox);
    label_2->setObjectName(QString::fromUtf8("label_2"));

    horizontalLayout_4->addWidget(label_2);

    b_intensity = new QDoubleSpinBox(groupBox);
    b_intensity->setObjectName(QString::fromUtf8("b_intensity"));
    b_intensity->setValue(1);

    horizontalLayout_4->addWidget(b_intensity);


    verticalLayout_3->addLayout(horizontalLayout_4);

    horizontalLayout_5 = new QHBoxLayout();
    horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
    label_3 = new QLabel(groupBox);
    label_3->setObjectName(QString::fromUtf8("label_3"));

    horizontalLayout_5->addWidget(label_3);

    c_intensity = new QDoubleSpinBox(groupBox);
    c_intensity->setObjectName(QString::fromUtf8("c_intensity"));
    c_intensity->setValue(1);

    horizontalLayout_5->addWidget(c_intensity);


    verticalLayout_3->addLayout(horizontalLayout_5);

    horizontalLayout_6 = new QHBoxLayout();
    horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
    label_4 = new QLabel(groupBox);
    label_4->setObjectName(QString::fromUtf8("label_4"));

    horizontalLayout_6->addWidget(label_4);

    x_intensity = new QDoubleSpinBox(groupBox);
    x_intensity->setObjectName(QString::fromUtf8("x_intensity"));
    x_intensity->setValue(1);

    horizontalLayout_6->addWidget(x_intensity);


    verticalLayout_3->addLayout(horizontalLayout_6);

    horizontalLayout_7 = new QHBoxLayout();
    horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
    label_5 = new QLabel(groupBox);
    label_5->setObjectName(QString::fromUtf8("label_5"));

    horizontalLayout_7->addWidget(label_5);

    y_intensity = new QDoubleSpinBox(groupBox);
    y_intensity->setObjectName(QString::fromUtf8("y_intensity"));
    y_intensity->setValue(1);

    horizontalLayout_7->addWidget(y_intensity);


    verticalLayout_3->addLayout(horizontalLayout_7);

    horizontalLayout_8 = new QHBoxLayout();
    horizontalLayout_8->setObjectName(QString::fromUtf8("horizontalLayout_8"));
    label_6 = new QLabel(groupBox);
    label_6->setObjectName(QString::fromUtf8("label_6"));

    horizontalLayout_8->addWidget(label_6);

    z_intensity = new QDoubleSpinBox(groupBox);
    z_intensity->setObjectName(QString::fromUtf8("z_intensity"));
    z_intensity->setValue(1);

    horizontalLayout_8->addWidget(z_intensity);


    verticalLayout_3->addLayout(horizontalLayout_8);

    horizontalLayout_9 = new QHBoxLayout();
    horizontalLayout_9->setObjectName(QString::fromUtf8("horizontalLayout_9"));
    label_7 = new QLabel(groupBox);
    label_7->setObjectName(QString::fromUtf8("label_7"));

    horizontalLayout_9->addWidget(label_7);

    rel_loss_intensity = new QSpinBox(groupBox);
    rel_loss_intensity->setObjectName(QString::fromUtf8("rel_loss_intensity"));
    rel_loss_intensity->setMaximum(100);
    rel_loss_intensity->setValue(10);

    horizontalLayout_9->addWidget(rel_loss_intensity);


    verticalLayout_3->addLayout(horizontalLayout_9);


    horizontalLayout_11->addWidget(groupBox);


    verticalLayout_5->addLayout(horizontalLayout_11);

    button_box = new QDialogButtonBox(TheoreticalSpectrumGenerationDialogTemplate);
    button_box->setObjectName(QString::fromUtf8("button_box"));
    QSizePolicy sizePolicy4(QSizePolicy::Ignored, QSizePolicy::Fixed);
    sizePolicy4.setHorizontalStretch(0);
    sizePolicy4.setVerticalStretch(0);
    sizePolicy4.setHeightForWidth(button_box->sizePolicy().hasHeightForWidth());
    button_box->setSizePolicy(sizePolicy4);
    button_box->setOrientation(Qt::Horizontal);
    button_box->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

    verticalLayout_5->addWidget(button_box);


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

#ifndef QT_NO_WHATSTHIS
    max_iso_label->setWhatsThis(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Defines the maximal isotopic peak which is added", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_WHATSTHIS

    max_iso_label->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Max Isotope:", 0, QApplication::UnicodeUTF8));
    generate_label->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Generate:", 0, QApplication::UnicodeUTF8));

    const bool __sortingEnabled = list_widget->isSortingEnabled();
    list_widget->setSortingEnabled(false);
    list_widget->item(0)->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "A-ions", 0, QApplication::UnicodeUTF8));
    list_widget->item(1)->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "B-ions", 0, QApplication::UnicodeUTF8));
    list_widget->item(2)->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "C-ions", 0, QApplication::UnicodeUTF8));
    list_widget->item(3)->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "X-ions", 0, QApplication::UnicodeUTF8));
    list_widget->item(4)->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Y-ions", 0, QApplication::UnicodeUTF8));
    list_widget->item(5)->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Z-ions", 0, QApplication::UnicodeUTF8));
    list_widget->item(6)->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Precursor", 0, QApplication::UnicodeUTF8));
    list_widget->item(7)->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Neutral losses", 0, QApplication::UnicodeUTF8));
    list_widget->item(8)->setText(QApplication::translate("TheoreticalSpectrumGenerationDialogTemplate", "Isotope clusters", 0, QApplication::UnicodeUTF8));

    list_widget->setSortingEnabled(__sortingEnabled);
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

QT_END_NAMESPACE

#endif // THEORETICALSPECTRUMGENERATIONDIALOGTEMPLATE_H

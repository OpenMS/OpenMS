// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------
//

#include <iostream>
#include <cmath>

#include <OpenMS/VISUAL/DIALOGS/SaveImageDialog.h>

#include <OpenMS/MATH/MISC/MathFunctions.h>

// Qt
#include <QtWidgets/QLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLabel>
#include <QValidator>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QImageWriter>
#include <QtWidgets/QApplication>

namespace OpenMS
{


  SaveImageDialog::SaveImageDialog(QWidget * parent) :
    QDialog(parent)
  {
    size_ratio_ = 1;
    //create dialog and layout (grid)
    QGridLayout * grid = new QGridLayout(this);

    //add accept/cancel buttons (and their layout)
    QBoxLayout * box_layout = new QHBoxLayout();
    grid->addLayout(box_layout, 5, 1);

    QPushButton * button = new QPushButton(this);
    button->setText("Cancel");
    connect(button, SIGNAL(clicked()), this, SLOT(reject()));
    box_layout->addWidget(button);

    button = new QPushButton(this);
    button->setText("Accept");
    button->setDefault(true);
    connect(button, SIGNAL(clicked()), this, SLOT(checkSize()));
    box_layout->addWidget(button);

    //add picture format selector
    QLabel * label = new QLabel("Picture format:", this);
    grid->addWidget(label, 0, 0);
    format_ = new QComboBox(this);
    QList<QByteArray> list = QImageWriter::supportedImageFormats();
    for (int i = 0; i < list.size(); ++i)
    {
      format_->insertItem(i, list.at(i));
    }
    grid->addWidget(format_, 0, 1, Qt::AlignLeft);
    //set format to PNG/JPEG if available
    int png = -1;
    int jpeg = -1;
    for (int i = 0; i < format_->count(); i++)
    {
      if (format_->itemText(i) == "PNG")
      {
        png = i;
      }
      if (format_->itemText(i) == "JPEG")
      {
        jpeg = i;
      }
    }
    if (png != -1)
    {
      format_->setCurrentIndex(png);
    }
    else if (jpeg != -1)
    {
      format_->setCurrentIndex(jpeg);
    }

    //add size boxes and label (and their layout)
    label = new QLabel("Size (WxH):", this);
    grid->addWidget(label, 1, 0);

    QValidator * v = new QIntValidator(1, 10000, this);
    box_layout = new QHBoxLayout();
    grid->addLayout(box_layout, 1, 1);
    size_x_ = new QLineEdit(this);
    size_x_->setValidator(v);
    connect(size_x_, SIGNAL(textChanged(const QString &)), this, SLOT(xSizeChanged(const QString &)));
    box_layout->addWidget(size_x_);
    label = new QLabel("x", this);
    box_layout->addWidget(label);
    size_y_ = new QLineEdit(this);
    size_y_->setValidator(v);
    connect(size_y_, SIGNAL(textChanged(const QString &)), this, SLOT(ySizeChanged(const QString &)));
    box_layout->addWidget(size_y_);
    label = new QLabel("pixel", this);
    box_layout->addWidget(label);

    size_proportions_ = new QCheckBox("keep proportions", this);
    size_proportions_->setChecked(true);
    connect(size_proportions_, SIGNAL(toggled(bool)), this, SLOT(proportionsActivated(bool)));
    grid->addWidget(size_proportions_, 2, 1);
  }

  void SaveImageDialog::setSize(int x, int y)
  {
    QString * temp = new QString();
    temp->setNum(x);
    size_x_->setText(*temp);
    temp->setNum(y);
    size_y_->setText(*temp);
    setSizeRatio_(float(x) / float(y));
  }

  void SaveImageDialog::setSizeRatio_(float r)
  {
    if (r == 0.0)
    {
      size_ratio_ = 1.0;
    }
    else
    {
      size_ratio_ = r;
    }
  }

  void SaveImageDialog::xSizeChanged(const QString & s)
  {
    if (size_proportions_->isChecked() && size_x_ == qApp->focusWidget())
    {
      QString * temp = new QString();
      temp->setNum((int)Math::round(s.toInt() / size_ratio_));
      size_y_->setText(*temp);
    }
  }

  void SaveImageDialog::ySizeChanged(const QString & s)
  {
    if (size_proportions_->isChecked() && size_y_ == qApp->focusWidget())
    {
      QString * temp = new QString();
      temp->setNum((int)Math::round(s.toInt() * size_ratio_));
      size_x_->setText(*temp);
    }
  }

  void SaveImageDialog::proportionsActivated(bool state)
  {
    if (state == true)
    {
      setSizeRatio_(QString(size_x_->text()).toFloat() / QString(size_y_->text()).toFloat());
    }
  }

  void SaveImageDialog::checkSize()
  {
    int x = size_x_->text().toInt();
    int y = size_y_->text().toInt();
    if (x > 0 && y > 0)
    {
      accept();
    }
  }

  int SaveImageDialog::getXSize()
  {
    return size_x_->text().toInt();
  }

  int SaveImageDialog::getYSize()
  {
    return size_y_->text().toInt();
  }

  QString SaveImageDialog::getFormat()
  {
    return format_->currentText();
  }

} //namespace

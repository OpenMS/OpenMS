// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/VISUAL/DIALOGS/HistogramDialog.h>

#include <QtWidgets/QPushButton>
#include <QtWidgets/QGridLayout>

using namespace std;

namespace OpenMS
{
  using namespace Math;

  HistogramDialog::HistogramDialog(const Histogram<> & distribution, QWidget * parent) :
    QDialog(parent)
  {
    setWindowTitle("Intensity Distribution");

    //layout
    QGridLayout * layout = new QGridLayout(this);
    layout->setRowStretch(0, 100);

    //ok
    QPushButton * ok_button_ = new QPushButton("&Apply Filter", this);
    ok_button_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    connect(ok_button_, SIGNAL(clicked()), this, SLOT(accept()));
    layout->addWidget(ok_button_, 1, 1);

    //cancel
    QPushButton * cancel_button_ = new QPushButton("&Cancel", this);
    cancel_button_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    connect(cancel_button_, SIGNAL(clicked()), this, SLOT(reject()));
    layout->addWidget(cancel_button_, 1, 2);

    //distribution
    mw_ = new HistogramWidget(distribution, this);
    mw_->showSplitters(true);
    layout->addWidget(mw_, 0, 0, 1, 3);

    //resize dialog
    adjustSize();
  }

  HistogramDialog::~HistogramDialog() = default;

  float HistogramDialog::getLeftSplitter()
  {
    return mw_->getLeftSplitter();
  }

  float HistogramDialog::getRightSplitter()
  {
    return mw_->getRightSplitter();
  }

  void HistogramDialog::setLeftSplitter(float position)
  {
    mw_->setLeftSplitter(position);
  }

  void HistogramDialog::setRightSplitter(float position)
  {
    mw_->setRightSplitter(position);
  }

  void HistogramDialog::setLegend(const String & legend)
  {
    mw_->setLegend(legend);
  }

  void HistogramDialog::setLogMode(bool log_mode)
  {
    mw_->setLogMode(log_mode);
  }

} //namespace OpenMS

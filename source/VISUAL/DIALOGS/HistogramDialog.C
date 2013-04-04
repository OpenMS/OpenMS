// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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


#include <OpenMS/VISUAL/DIALOGS/HistogramDialog.h>

#include <QtGui/QPushButton>
#include <QtGui/QGridLayout>

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

  HistogramDialog::~HistogramDialog()
  {
  }

  Real HistogramDialog::getLeftSplitter()
  {
    return mw_->getLeftSplitter();
  }

  Real HistogramDialog::getRightSplitter()
  {
    return mw_->getRightSplitter();
  }

  void HistogramDialog::setLeftSplitter(Real position)
  {
    mw_->setLeftSplitter(position);
  }

  void HistogramDialog::setRightSplitter(Real position)
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

// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

// OpenMS includes
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/CONCEPT/Types.h>

//qt includes
#include <QtGui/QPainter>
#include <QtGui/QColorDialog>
#include <QtGui/QPaintEvent>
#include <QtGui/QMouseEvent>

using namespace std;

namespace OpenMS
{

  ColorSelector::ColorSelector(QWidget * parent) :
    QWidget(parent),
    color_(255, 255, 255)
  {
    setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
  }

  QSize ColorSelector::sizeHint() const
  {
    return QSize(15, 15);
  }

  ColorSelector::~ColorSelector()
  {

  }

  void ColorSelector::paintEvent(QPaintEvent * /*e*/)
  {
    Int size = std::min(width(), height());
    QPainter painter(this);
    painter.setPen(QColor(0, 0, 0));
    painter.drawRect(0, 0, size - 1, size - 1);
    painter.setPen(QColor(255, 255, 255));
    painter.drawRect(1, 1, size - 3, size - 3);

    painter.fillRect(2, 2, size - 4, size - 4, color_);
  }

  void ColorSelector::mousePressEvent(QMouseEvent * e)
  {
    if (e->button() != Qt::LeftButton)
    {
      e->ignore();
      return;
    }
    QColor tmp = QColorDialog::getColor(color_, this);
    if (tmp.isValid())
    {
      color_ = tmp;
      repaint();
    }
  }

  const QColor & ColorSelector::getColor()
  {
    return color_;
  }

  void ColorSelector::setColor(const QColor & col)
  {
    color_ = col;
    repaint();
  }

} //namespace OpenMS

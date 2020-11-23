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
// $Maintainer: Jihyung Kim, Timo Sachsenberg  $
// $Authors: Jihyung Kim, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DVerticalLineItem.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>

#include <QtCore/QPoint>
#include <QtGui/QPainter>

using namespace std;

namespace OpenMS
{

  Annotation1DVerticalLineItem::Annotation1DVerticalLineItem(const double& x, const QColor& color, const QString& text) :
      Annotation1DItem(text),
      x_(x),
      color_(color)
  {
  }

  void Annotation1DVerticalLineItem::draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped)
  {
    //translate mz/intensity to pixel coordinates
    QPoint start_p, end_p;
    canvas->dataToWidget(x_, 0, start_p, flipped, true);
    canvas->dataToWidget(x_, canvas->getDataRange().maxY(), end_p, flipped, true);

    // draw line
    painter.drawLine(start_p, end_p);

    // compute bounding box on the specified painter
    // TODO: implement proper bounding box calculation
    // currently not needed as we don't support selection or moving
    bounding_box_ = QRectF(QPointF(start_p.x(), start_p.y()), QPointF(end_p.x(), end_p.y() ));    

    // TODO: draw according to proper bounding box to support switching axis and flipping
    // 5 pixel to x() was added to give some space between the line and the text
    painter.drawText(start_p.x() + 5, 20.0, text_);
  }

  void Annotation1DVerticalLineItem::move(const PointType & delta)
  {
    x_ += delta.getX();
  }

  void Annotation1DVerticalLineItem::setPosition(const double & x)
  {
    x_ = x;
  }

  const double & Annotation1DVerticalLineItem::getPosition() const
  {
    return x_;
  }

  void Annotation1DVerticalLineItem::ensureWithinDataRange(Plot1DCanvas* const)
  { // TODO: add code when needed
  }

} //Namespace

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
// $Authors: Johannes Junker, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

#include <QtGui/QPainter>
#include <QtCore/QPoint>

namespace OpenMS
{

  Annotation1DTextItem::Annotation1DTextItem(const PointType & position, const QString & text, int flags) :
    Annotation1DItem(text),
    position_(position),
    flags_(flags)
  {
  }

  Annotation1DTextItem::Annotation1DTextItem(const Annotation1DTextItem & rhs) :
    Annotation1DItem(rhs)
  {
    position_ = rhs.getPosition();
    flags_ = rhs.getFlags();
  }

  Annotation1DTextItem::~Annotation1DTextItem()
  {
  }

  void Annotation1DTextItem::draw(Spectrum1DCanvas * const canvas, QPainter & painter, bool flipped)
  {
    //translate mz/intensity to pixel coordinates
    QPoint pos;
    canvas->dataToWidget(position_.getX(), position_.getY(), pos, flipped, true);

    // compute bounding box of text_item on the specified painter
    bounding_box_ = painter.boundingRect(QRectF(pos, pos), flags_, text_);

    painter.drawText(bounding_box_, flags_, text_);
    if (selected_)
    {
      drawBoundingBox_(painter);
    }
  }

  void Annotation1DTextItem::move(const PointType & delta)
  {
    position_.setX(position_.getX() + delta.getX());
    position_.setY(position_.getY() + delta.getY());
  }

  void Annotation1DTextItem::setPosition(const Annotation1DTextItem::PointType & position)
  {
    position_ = position;
  }

  const Annotation1DTextItem::PointType & Annotation1DTextItem::getPosition() const
  {
    return position_;
  }

  void Annotation1DTextItem::setFlags(int flags)
  {
    flags_ = flags;
  }

  int Annotation1DTextItem::getFlags() const
  {
    return flags_;
  }

  void Annotation1DTextItem::ensureWithinDataRange(Spectrum1DCanvas * const canvas)
  {
    DRange<3> data_range = canvas->getDataRange();

    CoordinateType x_pos = position_.getX();
    CoordinateType y_pos = position_.getY() * canvas->getPercentageFactor();

    if (x_pos < data_range.minPosition()[0])
    {
      position_.setX(data_range.minPosition()[0]);
    }
    if (x_pos > data_range.maxPosition()[0])
    {
      position_.setX(data_range.maxPosition()[0]);
    }
    if (y_pos < data_range.minPosition()[1])
    {
      position_.setY(data_range.minPosition()[1] / canvas->getPercentageFactor());
    }
    if (y_pos > data_range.maxPosition()[1])
    {
      position_.setY(data_range.maxPosition()[1] / canvas->getPercentageFactor());
    }
  }

} // namespace OpenMS

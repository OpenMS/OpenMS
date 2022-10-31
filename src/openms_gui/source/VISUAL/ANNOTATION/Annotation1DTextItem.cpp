// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
#include <OpenMS/VISUAL/Plot1DCanvas.h>

namespace OpenMS
{
  namespace
  {
    Annotation1DTextItem p({0, 0}, "test");
  }

  void Annotation1DTextItem::ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index)
  {
    canvas->pushIntoDataRange(position_, layer_index);
  }

  void Annotation1DTextItem::draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped)
  {
    // translate units to pixel coordinates
    QPoint pos_text;
    canvas->dataToWidget(position_, pos_text, flipped);

    // compute bounding box of text_item on the specified painter
    bounding_box_ = painter.boundingRect(QRectF(pos_text, pos_text), flags_, text_);

    painter.drawText(bounding_box_, flags_, text_);
    if (selected_)
    {
      drawBoundingBox_(painter);
    }
  }

  void Annotation1DTextItem::move(const PointXYType delta, const Gravitator& /*gr*/, const DimMapper<2>& /*dim_mapper*/)
  {
    position_ += delta;
  }
} // namespace OpenMS

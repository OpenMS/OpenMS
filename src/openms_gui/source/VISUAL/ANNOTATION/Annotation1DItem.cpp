// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>

#include <QtGui/QInputDialog>
#include <QtGui/QPainter>

namespace OpenMS
{

  Annotation1DItem::Annotation1DItem(const QString & text) :
    bounding_box_(),
    selected_(true),
    text_(text)
  {
  }

  Annotation1DItem::Annotation1DItem(const Annotation1DItem & rhs)
  {
    bounding_box_ = rhs.boundingBox();
    selected_ = rhs.isSelected();
    text_ = rhs.getText();
  }

  Annotation1DItem::~Annotation1DItem()
  {
  }

  void Annotation1DItem::drawBoundingBox_(QPainter & painter)
  {
    // draw additional filled rectangles to highlight bounding box of selected distance_item
    painter.fillRect((int)(bounding_box_.topLeft().x()) - 3, (int)(bounding_box_.topLeft().y()) - 3, 3, 3, painter.pen().color());
    painter.fillRect((int)(bounding_box_.topRight().x()), (int)(bounding_box_.topRight().y()) - 3, 3, 3, painter.pen().color());
    painter.fillRect((int)(bounding_box_.bottomRight().x()), (int)(bounding_box_.bottomRight().y()), 3, 3, painter.pen().color());
    painter.fillRect((int)(bounding_box_.bottomLeft().x()) - 3, (int)(bounding_box_.bottomLeft().y()), 3, 3, painter.pen().color());
  }

  const QRectF & Annotation1DItem::boundingBox() const
  {
    return bounding_box_;
  }

  void Annotation1DItem::setSelected(bool selected)
  {
    selected_ = selected;
  }

  bool Annotation1DItem::isSelected() const
  {
    return selected_;
  }

  void Annotation1DItem::setText(const QString & text)
  {
    text_ = text;
  }

  const QString & Annotation1DItem::getText() const
  {
    return text_;
  }

  bool Annotation1DItem::editText()
  {
    bool ok;
    QString text = QInputDialog::getText(nullptr, "Edit text", "Enter text:", QLineEdit::Normal, this->getText(), &ok);
    if (ok && !text.isEmpty())
    {
      if (text == getText()) return false;
      this->setText(text);
      return true;
    }
    return false;
  }

} //Namespace

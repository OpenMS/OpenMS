// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>

#include <QtWidgets/QInputDialog>
#include <QPainter>

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

  Annotation1DItem::~Annotation1DItem() = default;

  void Annotation1DItem::drawBoundingBox_(QPainter & painter)
  {
    // draw additional filled rectangles to highlight bounding box of selected distance_item
    painter.fillRect((int)(bounding_box_.topLeft().x()) - 3, (int)(bounding_box_.topLeft().y()) - 3, 3, 3, painter.pen().color());
    painter.fillRect((int)(bounding_box_.topRight().x()), (int)(bounding_box_.topRight().y()) - 3, 3, 3, painter.pen().color());
    painter.fillRect((int)(bounding_box_.bottomRight().x()), (int)(bounding_box_.bottomRight().y()), 3, 3, painter.pen().color());
    painter.fillRect((int)(bounding_box_.bottomLeft().x()) - 3, (int)(bounding_box_.bottomLeft().y()), 3, 3, painter.pen().color());
  }

  const QRectF& Annotation1DItem::boundingBox() const
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
      if (text == getText())
      {
        return false;
      }
      this->setText(text);
      return true;
    }
    return false;
  }

} //Namespace

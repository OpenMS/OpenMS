// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/MultiGradientSelector.h>

//qt includes
#include <QPainter>
#include <QtWidgets/QColorDialog>
#include <QPixmap>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QPaintEvent>
#include <QContextMenuEvent>
#include <QtWidgets/QMenu>

using namespace std;

namespace OpenMS
{

  MultiGradientSelector::MultiGradientSelector(QWidget * parent) :
    QWidget(parent),
    gradient_(),
    margin_(5),
    gradient_area_width_(0),
    lever_area_height_(17),
    selected_(-1),
    selected_color_(Qt::white)
  {
    setMinimumSize(250, 45);
    setFocusPolicy(Qt::ClickFocus);
    setToolTip("Click the lever area to add new levers. Levers can be removed with the DEL key. "
               "Double click a lever to change its color. Levers can be dragged.<BR><BR>"
               "In the context menu you can select default gradients and change the interplation mode."
               );
  }

  MultiGradientSelector::~MultiGradientSelector() = default;

  const MultiGradient & MultiGradientSelector::gradient() const
  {
    return gradient_;
  }

  MultiGradient & MultiGradientSelector::gradient()
  {
    return gradient_;
  }

  void MultiGradientSelector::paintEvent(QPaintEvent * /* e */)
  {
    static QPixmap pixmap = QPixmap(size());
    pixmap.fill(palette().window().color());

    //calculate gradient area width
    if (gradient_area_width_ == 0)
    {
      gradient_area_width_ = width() - 2 * margin_ - 2;
    }

    QPainter painter(&pixmap);

    //gradient field outline
    painter.setPen(QColor(0, 0, 0));
    painter.drawRect(margin_, margin_, width() - 2 * margin_, height() - 2 * margin_ - lever_area_height_);

    //draw gradient
    for (Int i = 0; i <= gradient_area_width_; ++i)
    {
      painter.setPen(gradient_.interpolatedColorAt(i, 0, gradient_area_width_));
      painter.drawLine(margin_ + 1 + i, margin_ + 1, margin_ + 1 + i, height() - margin_ - lever_area_height_ - 1);
    }

    //levers
    painter.setPen(QColor(0, 0, 0));
    for (UInt i = 0; i < (UInt)gradient_.size(); ++i)
    {
      Int pos = Int(float(gradient_.position(i)) / 100.0 * gradient_area_width_ + margin_ + 1);
      painter.drawRect(pos - 4, height() - margin_ - lever_area_height_ + 5, 9, 9);
      painter.drawLine(pos - 4, height() - margin_ - lever_area_height_ + 5, pos, height() - margin_ - lever_area_height_);
      painter.drawLine(pos, height() - margin_ - lever_area_height_, pos + 4, height() - margin_ - lever_area_height_ + 5);
      painter.fillRect(pos - 3, height() - margin_ - lever_area_height_ + 6, 8, 8, gradient_.color(i));

      //selected lever
      if (Int(gradient_.position(i)) == selected_)
      {
        painter.fillRect(pos - 2, height() - margin_ - lever_area_height_ + 3, 6, 3, QColor(0, 0, 0));
        painter.fillRect(pos - 1, height() - margin_ - lever_area_height_ + 1, 4, 3, QColor(0, 0, 0));
      }
    }

    QPainter painter2(this);
    painter2.drawPixmap(0, 0, pixmap);
  }

  void MultiGradientSelector::mousePressEvent(QMouseEvent * e)
  {
    if (e->button() != Qt::LeftButton)
    {
      e->ignore();
      return;
    }

    left_button_pressed_ = true;

    // select lever
    // Starting with the rightmost lever, the first lever that overlaps the mouse pointer is selected (only lever with higher index can overlap lever with lower index).
    for (Int i = (Int)gradient_.size() - 1; i >= 0; --i)
    {
      Int pos = Int(float(gradient_.position(i)) / 100.0 * gradient_area_width_ + margin_ + 1);

      // mouse pointer over lever?
      if (e->x() >= pos - 3 && e->x() <= pos + 4 && e->y() >= height() - margin_ - lever_area_height_ + 8 && e->y() <= height() - margin_ - lever_area_height_ + 15)
      {
        selected_ = gradient_.position(i);
        selected_color_ = gradient_.color(i);
        repaint();
        return;
      }
    }

    //create new lever
    if (e->x() >= margin_ && e->x() <= width() - margin_ && e->y() >= height() - margin_ - lever_area_height_ && e->y() <= height() - margin_)
    {
      Int pos = Int(100 * (e->x() - margin_) / float(gradient_area_width_));
      gradient_.insert(pos, selected_color_);
      selected_ = pos;
      repaint();
    }
  }

  void MultiGradientSelector::mouseMoveEvent(QMouseEvent * e)
  {
    if (left_button_pressed_ && selected_ > 0 && selected_ < 100) // don't move first or last lever
    {
      //inside lever area
      if (e->x() >= margin_ && e->x() <= width() - margin_ && e->y() >= height() - margin_ - lever_area_height_ && e->y() <= height() - margin_)
      {
        Int pos = Int(100 * (e->x() - margin_) / float(gradient_area_width_));
        if (pos != selected_ && !gradient_.exists(pos))        // lever has been moved AND no other lever at the new position?
        {
          gradient_.remove(selected_);
          gradient_.insert(pos, selected_color_);
          selected_ = pos;
          repaint();
        }
      }
    }
  }

  void MultiGradientSelector::mouseReleaseEvent(QMouseEvent * e)
  {
    if (e->button() != Qt::LeftButton)
    {
      e->ignore();
      return;
    }
    left_button_pressed_ = false;
  }

  void MultiGradientSelector::mouseDoubleClickEvent(QMouseEvent * e)
  {
    for (UInt i = 0; i < (UInt)gradient_.size(); ++i)
    {
      Int pos = Int(float(gradient_.position(i)) / 100.0 * gradient_area_width_ + margin_ + 1);
      if (e->x() >= pos - 3 && e->x() <= pos + 4 && e->y() >= height() - margin_ - lever_area_height_ + 8 && e->y() <= height() - margin_ - lever_area_height_ + 15)
      {
        gradient_.insert(gradient_.position(i), QColorDialog::getColor(gradient_.color(i), this));
        if (Int(gradient_.position(i)) == selected_)
        {
          selected_color_ = gradient_.color(i);
        }
        return;
      }
    }
  }

  void MultiGradientSelector::keyPressEvent(QKeyEvent * e)
  {
    if (e->key() == Qt::Key_Delete && selected_ > 0 && selected_ < 100)     // don't remove first or last lever)
    {
      gradient_.remove(selected_);
      selected_ = -1;
      selected_color_ = Qt::white;
      repaint();
    }
    else
    {
      e->ignore();
    }
  }

  void MultiGradientSelector::stairsInterpolation(bool state)
  {
    if (state)
    {
      gradient_.setInterpolationMode(MultiGradient::IM_STAIRS);
    }
    else
    {
      gradient_.setInterpolationMode(MultiGradient::IM_LINEAR);
    }
  }

  void MultiGradientSelector::setInterpolationMode(MultiGradient::InterpolationMode mode)
  {
    gradient_.setInterpolationMode(mode);
  }

  MultiGradient::InterpolationMode MultiGradientSelector::getInterpolationMode() const
  {
    return gradient_.getInterpolationMode();
  }

  void MultiGradientSelector::contextMenuEvent(QContextMenuEvent * e)
  {
    QMenu main_menu(this);
    //Default gradient
    QMenu * defaults = main_menu.addMenu("Default gradients");
    defaults->addAction("grey - yellow - red - purple - blue - black");
    defaults->addAction("grey - black");
    defaults->addAction("yellow - red - purple - blue - black");
    defaults->addAction("orange - red - purple - blue - black");
    defaults->addAction("yellow - orange - red");
    defaults->addSeparator();
    defaults->addAction("black");
    defaults->addAction("white");
    defaults->addAction("red");
    defaults->addAction("green");
    defaults->addAction("blue");
    defaults->addAction("magenta");
    defaults->addAction("turquoise");
    defaults->addAction("yellow");

    //Interploate/Stairs
    QMenu * inter = main_menu.addMenu("Interpolation");
    QAction * current = inter->addAction("None");
    if (gradient_.getInterpolationMode() == MultiGradient::IM_STAIRS)
      current->setEnabled(false);
    current = inter->addAction("Linear");
    if (gradient_.getInterpolationMode() == MultiGradient::IM_LINEAR)
      current->setEnabled(false);

    //Execute
    QAction * result;
    if ((result = main_menu.exec(e->globalPos())))
    {
      if (result->text() == "grey - yellow - red - purple - blue - black")
      {
        gradient_ = MultiGradient::getDefaultGradientLinearIntensityMode();
      }
      if (result->text() == "grey - black")
      {
        gradient_ = MultiGradient::getDefaultGradientLogarithmicIntensityMode();
      }
      else if (result->text() == "yellow - red - purple - blue - black")
      {
        gradient_.fromString("Linear|0,#ffea00;6,#ff0000;14,#aa00ff;23,#5500ff;100,#000000");
      }
      else if (result->text() == "orange - red - purple - blue - black")
      {
        gradient_.fromString("Linear|0,#ffaa00;6,#ff0000;14,#aa00ff;23,#5500ff;100,#000000");
      }
      else if (result->text() == "yellow - orange - red")
      {
        gradient_.fromString("Linear|0,#ffea00;6,#ffaa00;100,#ff0000");
      }
      else if (result->text() == "black")
      {
        gradient_.fromString("Linear|0,#000000;100,#000000");
      }
      else if (result->text() == "white")
      {
        gradient_.fromString("Linear|0,#FFFFFF;100,#FFFFFF");
      }
      else if (result->text() == "red")
      {
        gradient_.fromString("Linear|0,#ff0000;100,#ff0000");
      }
      else if (result->text() == "green")
      {
        gradient_.fromString("Linear|0,#00ff00;100,#00ff00");
      }
      else if (result->text() == "blue")
      {
        gradient_.fromString("Linear|0,#0000ff;100,#0000ff");
      }
      else if (result->text() == "magenta")
      {
        gradient_.fromString("Linear|0,#ff00ff;100,#ff00ff");
      }
      else if (result->text() == "turquoise")
      {
        gradient_.fromString("Linear|0,#00ffff;100,#00ffff");
      }
      else if (result->text() == "yellow")
      {
        gradient_.fromString("Linear|0,#ffff00;100,#ffff00");
      }
      else if (result->text() == "None")
      {
        setInterpolationMode(MultiGradient::IM_STAIRS);
      }
      else if (result->text() == "Linear")
      {
        setInterpolationMode(MultiGradient::IM_LINEAR);
      }
    }
  }

} //namespace OpenMS

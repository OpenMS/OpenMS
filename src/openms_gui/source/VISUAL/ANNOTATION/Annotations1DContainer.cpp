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

#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>

#include <QtCore/QPoint>
#include <QtCore/QObject>
#include <QtCore/QRectF>
#include <QtGui/QPainter>

#include <iostream>

namespace OpenMS
{

  Annotations1DContainer::Annotations1DContainer() :
    std::list<Annotation1DItem *>()
  {
  }

  Annotations1DContainer::Annotations1DContainer(const Annotations1DContainer & rhs) :
    std::list<Annotation1DItem *>()
  {
    //copy annotations
    Annotation1DItem * new_item = nullptr;
    for (ConstIterator it = rhs.begin(); it != rhs.end(); ++it)
    {
      const Annotation1DDistanceItem * distance_item = dynamic_cast<const Annotation1DDistanceItem *>(*it);
      if (distance_item)
      {
        new_item = new Annotation1DDistanceItem(*distance_item);
        push_back(new_item);
        continue;
      }
      const Annotation1DTextItem * text_item = dynamic_cast<const Annotation1DTextItem *>(*it);
      if (text_item)
      {
        new_item = new Annotation1DTextItem(*text_item);
        push_back(new_item);
        continue;
      }
      const Annotation1DPeakItem * peak_item = dynamic_cast<const Annotation1DPeakItem *>(*it);
      if (peak_item)
      {
        new_item = new Annotation1DPeakItem(*peak_item);
        push_back(new_item);
        continue;
      }
    }
  }

  Annotations1DContainer & Annotations1DContainer::operator=(const Annotations1DContainer & rhs)
  {
    if (this != &rhs)
    {
      //delete existing annotations
      for (Iterator it = begin(); it != end(); ++it)
      {
        delete *it;
      }
      //clear list
      clear();
      //copy annotations
      Annotation1DItem * new_item = nullptr;
      for (ConstIterator it = rhs.begin(); it != rhs.end(); ++it)
      {
        const Annotation1DDistanceItem * distance_item = dynamic_cast<const Annotation1DDistanceItem *>(*it);
        if (distance_item)
        {
          new_item = new Annotation1DDistanceItem(*distance_item);
          push_back(new_item);
          continue;
        }
        const Annotation1DTextItem * text_item = dynamic_cast<const Annotation1DTextItem *>(*it);
        if (text_item)
        {
          new_item = new Annotation1DTextItem(*text_item);
          push_back(new_item);
          continue;
        }
        const Annotation1DPeakItem * peak_item = dynamic_cast<const Annotation1DPeakItem *>(*it);
        if (peak_item)
        {
          new_item = new Annotation1DPeakItem(*peak_item);
          push_back(new_item);
          continue;
        }
      }
    }
    return *this;
  }

  Annotations1DContainer::~Annotations1DContainer()
  {
    for (Iterator it = begin(); it != end(); ++it)
    {
      delete *it;
    }
  }

  Annotation1DItem * Annotations1DContainer::getItemAt(const QPoint & pos) const
  {
    for (ConstIterator it = begin(); it != end(); ++it)
    {
      if ((*it)->boundingBox().contains(pos))
      {
        return *it;
      }
    }
    return nullptr;
  }

  void Annotations1DContainer::selectItemAt(const QPoint & pos)
  {
    Annotation1DItem * item = getItemAt(pos);
    if (item != nullptr)
    {
      item->setSelected(true);
    }
  }

  void Annotations1DContainer::deselectItemAt(const QPoint & pos)
  {
    Annotation1DItem * item = getItemAt(pos);
    if (item != nullptr)
    {
      item->setSelected(false);
    }
  }

  void Annotations1DContainer::selectAll()
  {
    for (Iterator it = begin(); it != end(); ++it)
    {
      (*it)->setSelected(true);
    }
  }

  void Annotations1DContainer::deselectAll()
  {
    for (Iterator it = begin(); it != end(); ++it)
    {
      (*it)->setSelected(false);
    }
  }

  void Annotations1DContainer::removeSelectedItems()
  {
    for (Iterator it = begin(); it != end(); )
    {
      if ((*it)->isSelected())
      {
        delete *it;
        it = erase(it);
      }
      else
      {
        ++it;
      }
    }
  }

  std::vector<Annotation1DItem*> Annotations1DContainer::getSelectedItems()
  {
    // initialize with maximal possible size
    std::vector<Annotation1DItem*> annotation_items(size());
    // copy if is selected
    auto it = std::copy_if(begin(), end(), annotation_items.begin(), [](Annotation1DItem* anno){return anno->isSelected();});
    // resize to number of actually copied items
    annotation_items.resize(std::distance(annotation_items.begin(), it));

    return annotation_items;
  }

} //Namespace

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>


namespace OpenMS
{
  Annotations1DContainer::Annotations1DContainer() = default;

  Annotations1DContainer::Annotations1DContainer(const Annotations1DContainer& rhs)
    : Base() // do not copy items here. They need cloning!
  {
    // copy annotations
    for (auto item : rhs)
    {
      push_back(item->clone());
    }
  }

  Annotations1DContainer& Annotations1DContainer::operator=(const Annotations1DContainer & rhs)
  {
    if (this != &rhs)
    {
      deleteAllItems_();
      // clear list
      clear();

      // copy annotations
      for (auto ptr_item : rhs)
      {
        push_back(ptr_item->clone());
      }
    }
    return *this;
  }


  Annotations1DContainer::~Annotations1DContainer()
  {
    deleteAllItems_();
  }

  Annotation1DItem * Annotations1DContainer::getItemAt(const QPoint & pos) const
  {
    for (const auto ptr_item : *this)
    {
      if (ptr_item->boundingBox().contains(pos))
      {
        return ptr_item;
      }
    }
    return nullptr;
  }

  void Annotations1DContainer::selectItemAt(const QPoint & pos) const
  {
    Annotation1DItem* item = getItemAt(pos);
    if (item != nullptr)
    {
      item->setSelected(true);
    }
  }

  void Annotations1DContainer::deselectItemAt(const QPoint & pos) const
  {
    Annotation1DItem* item = getItemAt(pos);
    if (item != nullptr)
    {
      item->setSelected(false);
    }
  }

  void Annotations1DContainer::selectAll()
  {
    for (const auto ptr_item : *this)
    {
      ptr_item->setSelected(true);
    }
  }

  void Annotations1DContainer::deselectAll()
  {
    for (const auto ptr_item : *this)
    {
      ptr_item->setSelected(false);
    }
  }

  void Annotations1DContainer::removeSelectedItems()
  {
    for (auto it = begin(); it != end();)
    {
      if ((*it)->isSelected())
      {
        delete (*it);
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

  
  void Annotations1DContainer::deleteAllItems_() const
  {
    for (const auto ptr_item : *this)
    {
      delete ptr_item;
    }
  }

} //Namespace

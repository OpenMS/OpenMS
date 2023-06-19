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

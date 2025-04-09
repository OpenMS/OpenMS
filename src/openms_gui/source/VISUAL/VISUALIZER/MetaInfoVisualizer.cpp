// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoVisualizer.h>

//QT
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QButtonGroup>

//STL
#include <iostream>
#include <vector>
#include <utility>

using namespace std;

namespace OpenMS
{

  MetaInfoVisualizer::MetaInfoVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<MetaInfoInterface>()
  {
    buttongroup_ = new QButtonGroup();
    nextrow_ = 0;
    viewlayout_ = new QGridLayout();

    addLabel_("Modify MetaData information.");
    addSeparator_();
    mainlayout_->addLayout(viewlayout_, row_, 0, 1, 3);
    //increase row counter for mainlayout_.
    row_++;
  }

  void MetaInfoVisualizer::load(MetaInfoInterface & m)
  {
    ptr_ = &m;
    temp_ = m;

    //keys_ is a vector of indices of type UInt
    temp_.getKeys(keys_);

    //Load actual metaInfo Data into viewLayout_
    for (Size i = 0; i < keys_.size(); ++i)
    {
      loadData_(keys_[i]);
    }

    addSeparator_();
    addLabel_("Add new MetaInfo entry.");
    addLineEdit_(newkey_, "Key");
    addLineEdit_(newdescription_, "Description");
    addLineEdit_(newvalue_, "Value");
    add2Buttons_(addbutton_, "Add", clearbutton_, "Clear");
    if (!isEditable())
    {
      addbutton_->setEnabled(false);
      clearbutton_->setEnabled(false);
    }

    finishAdding_();

    connect(buttongroup_, SIGNAL(buttonClicked(int)), this, SLOT(remove_(int)));
    connect(addbutton_, SIGNAL(clicked()), this, SLOT(add_()));
    connect(clearbutton_, SIGNAL(clicked()), this, SLOT(clear_()));
  }

  //----------------------------------------------------------------------------
  // SLOTS
  //----------------------------------------------------------------------------

  void MetaInfoVisualizer::remove_(int index)
  {
    UInt id = (UInt)index;

    //Remove label
    std::vector<std::pair<UInt, QLabel *> >::iterator iter;
    for (iter = metalabels_.begin(); iter < metalabels_.end(); ++iter)
    {
      if ((*iter).first == id)
      {
        viewlayout_->removeWidget((*iter).second);
        (*iter).second->hide();
        (*iter).second = 0;
        delete(*iter).second;
        // cppcheck-suppress eraseDereference
        metalabels_.erase(iter);
      }
    }

    //Remove QLineEdit
    std::vector<std::pair<UInt, QLineEdit *> >::iterator iter2;
    for (iter2 = metainfoptr_.begin(); iter2 < metainfoptr_.end(); ++iter2)
    {
      if ((*iter2).first == id)
      {
        viewlayout_->removeWidget((*iter2).second);
        (*iter2).second->hide();
        (*iter2).second = 0;
        delete(*iter2).second;
        // cppcheck-suppress eraseDereference
        metainfoptr_.erase(iter2);
      }


    }

    //Remove QButton
    std::vector<std::pair<UInt, QAbstractButton *> >::iterator iter3 = metabuttons_.begin();
    while (iter3 != metabuttons_.end())
    {
      if ((*iter3).first == id)
      {
        viewlayout_->removeWidget((*iter3).second);
        (*iter3).second->hide();
        (*iter3).second = 0;
        delete(*iter3).second;
        iter3 = metabuttons_.erase(iter3);
      }
      else
      {
        ++iter3;
      }
    }

    //Remove entry from metainfointerface
    temp_.removeMetaValue(id);
    temp_.getKeys(keys_);
  }

  void MetaInfoVisualizer::loadData_(UInt index)
  {
    //----------------------------------------------------------------------------
    //  All metainfo goes into the viewlayout_
    //----------------------------------------------------------------------------

    QLabel * lab;
    QLineEdit * ptr;
    QPushButton * button;


    lab = new QLabel(temp_.metaRegistry().getName(index).c_str(), this);
    viewlayout_->addWidget(lab, nextrow_, 0);

    ptr = new QLineEdit(this);
    ptr->setText(temp_.getMetaValue(index).toString().c_str());
    viewlayout_->addWidget(ptr, nextrow_, 1);

    button = new QPushButton("Remove", this);
    if (!isEditable())
    {
      button->setEnabled(false);
    }
    viewlayout_->addWidget(button, nextrow_, 2);

    //Store information about ID(index) and QWidget
    metalabels_.emplace_back(index, lab);
    metainfoptr_.emplace_back(index, ptr);
    metabuttons_.emplace_back(index, button);

    //Insert new button with ID into buttongroup
    buttongroup_->addButton(button, index);
    nextrow_++;

    lab->show();
    ptr->show();
    button->show();
  }

  void MetaInfoVisualizer::add_()
  {
    String name(newkey_->text());
    String description(newdescription_->text());
    String value(newvalue_->text());


    if (name.trim().length() == 0)    //Must have a name
    {
      return;
    }

    //Register new entry and update metainfointerface object
    UInt newindex = temp_.metaRegistry().registerName(name, description, "");

    //Store new data in temporary metainfo object
    temp_.setMetaValue(newindex, value);
    temp_.getKeys(keys_);


    //Update viewlayout_
    try
    {
      //check whether there is already an entry in GUI for added metainfo.
      //If index already exists, return and do nothing.
      if (buttongroup_->button(newindex) != nullptr)
      {
        return;
      }

      loadData_(newindex);

    }
    catch (exception & e)
    {
      std::cout << "Error while trying to create new meta info data. " << e.what() << endl;
    }
  }

  void MetaInfoVisualizer::clear_()
  {
    newkey_->clear();
    newdescription_->clear();
    newvalue_->clear();
  }

  void MetaInfoVisualizer::store()
  {
    //Store QLineEdit information
    std::vector<std::pair<UInt, QLineEdit *> >::iterator iter2;
    for (iter2 = metainfoptr_.begin(); iter2 < metainfoptr_.end(); ++iter2)
    {
      UInt index = (*iter2).first;
      String value(((*iter2).second)->text());
      temp_.setMetaValue(index, value);
    }
    //copy temporary stored data into metainfo object
    (*ptr_) = temp_;

  }

  void MetaInfoVisualizer::undo_()
  {
    try
    {
      //Delete all data in GUI
      std::vector<UInt> keys_temp = keys_;
      for (Size i = 0; i < keys_temp.size(); ++i)
      {
        remove_(keys_temp[i]);
      }

      metalabels_.clear();
      metainfoptr_.clear();
      metabuttons_.clear();

      //Restore original data and refill GUI
      temp_ = (*ptr_);
      nextrow_ = 0;
      keys_.clear();
      ptr_->getKeys(keys_);
      for (Size i = 0; i < keys_.size(); ++i)
      {
        loadData_(keys_[i]);
      }
    }
    catch (exception & e)
    {
      cout << "Error while trying to restore original metainfo data. " << e.what() << endl;
    }
  }

}

// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  this library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/SoftwareVisualizer.h>

//QT
#include <QtGui/QTextEdit>
#include <QtGui/QLineEdit>

#include <iostream>

using namespace std;

namespace OpenMS
{

  SoftwareVisualizer::SoftwareVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Software>()
  {
    addLabel_("Modify software information.");
    addSeparator_();
    addLineEdit_(software_name_, "Name");
    addLineEdit_(software_version_, "Version");

    finishAdding_();
  }

  void SoftwareVisualizer::update_()
  {
    software_name_->setText(temp_.getName().c_str());
    software_version_->setText(temp_.getVersion().c_str());
  }

  void SoftwareVisualizer::store()
  {
    ptr_->setName(software_name_->text());
    ptr_->setVersion(software_version_->text());

    temp_ = (*ptr_);
  }

  void SoftwareVisualizer::undo_()
  {
    update_();
  }

}

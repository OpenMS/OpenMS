// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <OpenMS/DATASTRUCTURES/Param.h>


// QT includes
#include <QtGui/QButtonGroup>

// STL includes
#include <iostream>

using namespace std;

namespace OpenMS
{

  TOPPViewOpenDialog::TOPPViewOpenDialog(const String & data_name, bool as_window, bool as_2d, bool cutoff, QWidget * parent) :
    QDialog(parent),
    map_as_2d_disabled_(false)
  {
    setupUi(this);

    //init map view
    QButtonGroup * button_group = new QButtonGroup(this);
    button_group->addButton(d1_);
    button_group->addButton(d2_);
    button_group->addButton(d3_);
    if (!as_2d)
    {
      d1_->setChecked(true);
      d1_->setFocus();
    }
    else
    {
      d2_->setChecked(true);
      d2_->setFocus();
    }

    //init intensity cutoff
    button_group = new QButtonGroup(this);
    button_group->addButton(cutoff_);
    button_group->addButton(nocutoff_);
    if (!cutoff)
    {
      nocutoff_->setChecked(true);
      cutoff_->setFocus();
    }
    else
    {
      cutoff_->setChecked(true);
      cutoff_->setFocus();
    }

    //init open as
    button_group = new QButtonGroup(this);
    button_group->addButton(window_);
    button_group->addButton(layer_);
    button_group->addButton(merge_);
    connect(button_group, SIGNAL(buttonClicked(QAbstractButton *)), this, SLOT(updateViewMode_(QAbstractButton *)));
    if (!as_window)
    {
      layer_->setChecked(true);
      layer_->setFocus();
    }
    else
    {
      window_->setChecked(true);
      window_->setFocus();
    }
    connect(merge_combo_, SIGNAL(activated(int)), merge_, SLOT(click()));

    //set title
    setWindowTitle((String("Open data options for ") + data_name).toQString());
  }

  TOPPViewOpenDialog::~TOPPViewOpenDialog()
  {
  }

  bool TOPPViewOpenDialog::viewMapAs2D() const
  {
    if (d2_->isChecked())
      return true;

    return false;
  }

  bool TOPPViewOpenDialog::viewMapAs1D() const
  {
    if (d1_->isChecked())
      return true;

    return false;
  }

  bool TOPPViewOpenDialog::isCutoffEnabled() const
  {
    if (cutoff_->isChecked())
      return true;

    return false;
  }

  bool TOPPViewOpenDialog::openAsNewWindow() const
  {
    if (window_->isChecked())
      return true;

    return false;
  }

  void TOPPViewOpenDialog::disableDimension(bool as_2d)
  {
    d1_->setChecked(!as_2d);
    d1_->setEnabled(false);
    d2_->setChecked(as_2d);
    d2_->setEnabled(false);
    d3_->setEnabled(false);
    map_as_2d_disabled_ = true;
  }

  void TOPPViewOpenDialog::disableCutoff(bool cutoff_on)
  {
    cutoff_->setChecked(cutoff_on);
    cutoff_->setEnabled(false);
    nocutoff_->setEnabled(false);
  }

  void TOPPViewOpenDialog::disableLocation(bool as_window)
  {
    window_->setEnabled(false);
    layer_->setEnabled(false);
    merge_->setEnabled(false);
    merge_combo_->setEnabled(false);
    if (as_window)
    {
      window_->setChecked(true);
    }
    else
    {
      layer_->setChecked(true);
    }
  }

  void TOPPViewOpenDialog::updateViewMode_(QAbstractButton * button)
  {
    if (button == layer_ || button == merge_)
    {
      d1_->setEnabled(false);
      d2_->setEnabled(false);
      d3_->setEnabled(false);
    }
    else if (!map_as_2d_disabled_)
    {
      d1_->setEnabled(true);
      d2_->setEnabled(true);
      d3_->setEnabled(true);
    }
  }

  void TOPPViewOpenDialog::setMergeLayers(const Map<Size, String> & layers)
  {
    //remove all items
    merge_combo_->clear();

    if (layers.size() != 0)
    {
      merge_->setEnabled(true);
      merge_combo_->setEnabled(true);
      UInt i = 0;
      for (Map<Size, String>::const_iterator it = layers.begin(); it != layers.end(); ++it)
      {
        merge_combo_->insertItem(i++, it->second.toQString(), (int)(it->first));
      }
    }
    else
    {
      merge_->setEnabled(false);
      merge_combo_->setEnabled(false);
    }
  }

  Int TOPPViewOpenDialog::getMergeLayer() const
  {
    if (merge_->isChecked())
    {
      return merge_combo_->itemData(merge_combo_->currentIndex()).toInt();
    }

    return -1;
  }

}

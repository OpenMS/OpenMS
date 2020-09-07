// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <ui_TOPPViewOpenDialog.h>

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/Param.h>


// QT includes
#include <QtWidgets/QButtonGroup>

// STL includes
#include <iostream>

using namespace std;


namespace OpenMS
{

  TOPPViewOpenDialog::TOPPViewOpenDialog(const String & data_name, bool as_window, bool as_2d, bool cutoff, QWidget * parent) :
    QDialog(parent),
    map_as_2d_disabled_(false),
    ui_(new Ui::TOPPViewOpenDialogTemplate)
  {
    ui_->setupUi(this);

    //init map view
    if (!as_2d)
    {
      ui_->d1_->setChecked(true);
      ui_->d1_->setFocus();
    }
    else
    {
      ui_->d2_->setChecked(true);
      ui_->d2_->setFocus();
    }

    //init intensity cutoff
    if (cutoff)
    {
      ui_->intensity_cutoff_->setChecked(true);
    }

    //init open as
    if (!as_window)
    {
      ui_->layer_->setChecked(true);
      ui_->layer_->setFocus();
    }
    else
    {
      ui_->window_->setChecked(true);
      ui_->window_->setFocus();
    }
    connect(ui_->merge_combo_, SIGNAL(activated(int)), ui_->merge_, SLOT(click()));

    //set title
    setWindowTitle((String("Open data options for ") + data_name).toQString());
  }

  TOPPViewOpenDialog::~TOPPViewOpenDialog()
  {
    delete ui_;
  }

  bool TOPPViewOpenDialog::viewMapAs2D() const
  {
    return ui_->d2_->isChecked();
  }

  bool TOPPViewOpenDialog::viewMapAs1D() const
  {
    return ui_->d1_->isChecked();
  }

  bool TOPPViewOpenDialog::isCutoffEnabled() const
  {
    return ui_->intensity_cutoff_->isChecked();
  }

  bool TOPPViewOpenDialog::isDataDIA() const
  {
    return ui_->dia_data_->isChecked();
  }

  bool TOPPViewOpenDialog::openAsNewWindow() const
  {
    return ui_->window_->isChecked();
  }

  void TOPPViewOpenDialog::disableDimension(bool as_2d)
  {
    ui_->d1_->setChecked(!as_2d);
    ui_->d1_->setEnabled(false);
    ui_->d2_->setChecked(as_2d);
    ui_->d2_->setEnabled(false);
    ui_->d3_->setEnabled(false);
    map_as_2d_disabled_ = true;
  }

  void TOPPViewOpenDialog::disableCutoff(bool /* cutoff_on */)
  {
    ui_->intensity_cutoff_->setChecked(false);
  }

  void TOPPViewOpenDialog::disableLocation(bool as_window)
  {
    ui_->window_->setEnabled(false);
    ui_->layer_->setEnabled(false);
    ui_->merge_->setEnabled(false);
    ui_->merge_combo_->setEnabled(false);
    if (as_window)
    {
      ui_->window_->setChecked(true);
    }
    else
    {
      ui_->layer_->setChecked(true);
    }
  }

  void TOPPViewOpenDialog::updateViewMode_(QAbstractButton * button)
  {
    if (button == ui_->layer_ || button == ui_->merge_)
    {
      ui_->d1_->setEnabled(false);
      ui_->d2_->setEnabled(false);
      ui_->d3_->setEnabled(false);
    }
    else if (!map_as_2d_disabled_)
    {
      ui_->d1_->setEnabled(true);
      ui_->d2_->setEnabled(true);
      ui_->d3_->setEnabled(true);
    }
  }

  void TOPPViewOpenDialog::setMergeLayers(const Map<Size, String> & layers)
  {
    // remove all items
    ui_->merge_combo_->clear();

    if (layers.size() != 0)
    {
      ui_->merge_->setEnabled(true);
      ui_->merge_combo_->setEnabled(true);
      UInt i = 0;
      for (Map<Size, String>::const_iterator it = layers.begin(); it != layers.end(); ++it)
      {
        ui_->merge_combo_->insertItem(i++, it->second.toQString(), (int)(it->first));
      }
    }
    else
    {
      ui_->merge_->setEnabled(false);
      ui_->merge_combo_->setEnabled(false);
    }
  }

  Int TOPPViewOpenDialog::getMergeLayer() const
  {
    if (ui_->merge_->isChecked())
    {
      return ui_->merge_combo_->itemData(ui_->merge_combo_->currentIndex()).toInt();
    }

    return -1;
  }

}

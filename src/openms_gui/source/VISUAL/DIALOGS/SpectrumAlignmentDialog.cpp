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
#include <OpenMS/VISUAL/DIALOGS/SpectrumAlignmentDialog.h>
#include <ui_SpectrumAlignmentDialog.h>


#include <OpenMS/VISUAL/Spectrum1DWidget.h>

// QT includes
#include <QtWidgets/QButtonGroup>

#include <vector>

namespace OpenMS
{
  SpectrumAlignmentDialog::SpectrumAlignmentDialog(Spectrum1DWidget * parent) :
    layer_indices_1_(),
    layer_indices_2_(),
    ui_(new Ui::SpectrumAlignmentDialogTemplate)
  {
    ui_->setupUi(this);

    QButtonGroup * button_group = new QButtonGroup(this);
    button_group->addButton(ui_->ppm);
    button_group->addButton(ui_->da);
    ui_->da->setChecked(true);

    Spectrum1DCanvas * cc = parent->canvas();
    for (UInt i = 0; i < cc->getLayerCount(); ++i)
    {
      const LayerData & layer = cc->getLayer(i);
      if (layer.flipped)
      {
        ui_->layer_list_2->addItem(layer.getName().toQString());
        layer_indices_2_.push_back(i);
      }
      else
      {
        ui_->layer_list_1->addItem(layer.getName().toQString());
        layer_indices_1_.push_back(i);
      }
    }
    // select first item of each list
    if (ui_->layer_list_1->count() > 0)
    {
      ui_->layer_list_1->setCurrentRow(0);
    }
    if (ui_->layer_list_2->count() > 0)
    {
      ui_->layer_list_2->setCurrentRow(0);
    }
  }

  SpectrumAlignmentDialog::~SpectrumAlignmentDialog()
  {
    delete ui_;
  }

  double SpectrumAlignmentDialog::getTolerance() const
  {
    return ui_->tolerance_spinbox->value();
  }

  bool SpectrumAlignmentDialog::isPPM() const
  {
    return ui_->ppm->isChecked();
  }


  Int SpectrumAlignmentDialog::get1stLayerIndex()
  {
    if (ui_->layer_list_1->count() == 0 || ui_->layer_list_1->currentRow() == -1)
    {
      return -1;
    }
    if (layer_indices_1_.size() > (Size)(ui_->layer_list_1->currentRow()))
    {
      return layer_indices_1_[(Size)(ui_->layer_list_1->currentRow())];
    }
    return -1;
  }

  Int SpectrumAlignmentDialog::get2ndLayerIndex()
  {
    if (ui_->layer_list_2->count() == 0 || ui_->layer_list_2->currentRow() == -1)
    {
      return -1;
    }
    if (layer_indices_2_.size() > (Size)(ui_->layer_list_2->currentRow()))
    {
      return layer_indices_2_[(Size)(ui_->layer_list_2->currentRow())];
    }
    return -1;
  }

} // namespace

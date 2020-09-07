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
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>
#include <ui_Spectrum2DGoToDialog.h>


#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtWidgets/QLineEdit>

using namespace std;

namespace OpenMS
{

  Spectrum2DGoToDialog::Spectrum2DGoToDialog(QWidget * parent) :
    QDialog(parent),
    ui_(new Ui::Spectrum2DGoToDialogTemplate)
  {
    ui_->setupUi(this);
  }

  Spectrum2DGoToDialog::~Spectrum2DGoToDialog()
  {
    delete ui_;
  }

  void Spectrum2DGoToDialog::setRange(float min_rt, float max_rt, float min_mz, float max_mz)
  {
    ui_->min_rt_->setText(QString::number(min_rt));
    ui_->max_rt_->setText(QString::number(max_rt));
    ui_->min_mz_->setText(QString::number(min_mz));
    ui_->max_mz_->setText(QString::number(max_mz));
  }

  void Spectrum2DGoToDialog::setMinMaxOfRange(float min_rt, float max_rt, float min_mz, float max_mz)
  {
    ui_->min_rt_const_->setText(QString("min: ") + QString::number(min_rt));
    ui_->max_rt_const_->setText(QString("max: ") + QString::number(max_rt));
    ui_->min_mz_const_->setText(QString("min: ") + QString::number(min_mz));
    ui_->max_mz_const_->setText(QString("max: ") + QString::number(max_mz));
  }

  bool Spectrum2DGoToDialog::checked()
  {
    return ui_->clip_checkbox->checkState() == Qt::Checked;
  }

  void Spectrum2DGoToDialog::fixRange()
  {
    // load from GUI
    float min_rt = ui_->min_rt_->text().toFloat();
    float max_rt = ui_->max_rt_->text().toFloat();
    float min_mz = ui_->min_mz_->text().toFloat();
    float max_mz = ui_->max_mz_->text().toFloat();

    // ensure correct order of min and max
    if (min_rt > max_rt) swap(min_rt, max_rt);
    if (min_mz > max_mz) swap(min_mz, max_mz);

    // do not allow range of 0 --> extend to 1 sec
    if (min_rt == max_rt)
    {
      min_rt -= 0.5;
      max_rt += 0.5;
    }
    if (min_mz == max_mz)
    {
      min_mz -= 0.5;
      max_mz += 0.5;
    }

    // store in GUI
    ui_->min_rt_->setText(QString::number(min_rt));
    ui_->max_rt_->setText(QString::number(max_rt));
    ui_->min_mz_->setText(QString::number(min_mz));
    ui_->max_mz_->setText(QString::number(max_mz));
  }

  float Spectrum2DGoToDialog::getMinRT() const
  {
    return ui_->min_rt_->text().toFloat();
  }

  float Spectrum2DGoToDialog::getMaxRT() const
  {
    return ui_->max_rt_->text().toFloat();
  }

  float Spectrum2DGoToDialog::getMinMZ() const
  {
    return ui_->min_mz_->text().toFloat();
  }

  float Spectrum2DGoToDialog::getMaxMZ() const
  {
    return ui_->max_mz_->text().toFloat();
  }

  void Spectrum2DGoToDialog::enableFeatureNumber(bool enabled)
  {
    ui_->feature_label_->setEnabled(enabled);
    ui_->nr_->setEnabled(enabled);
    ui_->feature_number_->setEnabled(enabled);
    //Reorder tab order
    if (enabled)
    {
      setTabOrder(ui_->feature_number_, ui_->ok_button_);
      setTabOrder(ui_->ok_button_, ui_->cancel_button_);
      setTabOrder(ui_->cancel_button_, ui_->min_mz_);
      setTabOrder(ui_->min_mz_, ui_->max_mz_);
      setTabOrder(ui_->max_mz_, ui_->min_rt_);
      setTabOrder(ui_->min_rt_, ui_->max_rt_);
    }
    else
    {
      setTabOrder(ui_->min_mz_, ui_->max_mz_);
      setTabOrder(ui_->max_mz_, ui_->min_rt_);
      setTabOrder(ui_->min_rt_, ui_->max_rt_);
      setTabOrder(ui_->max_rt_, ui_->ok_button_);
      setTabOrder(ui_->ok_button_, ui_->cancel_button_);
    }
  }

  String Spectrum2DGoToDialog::getFeatureNumber() const
  {
    return ui_->feature_number_->text();
  }

  bool Spectrum2DGoToDialog::showRange() const
  {
    if (ui_->feature_number_->text().trimmed() != "")
    {
      return false;
    }
    return true;
  }

} //namespace OpenMS

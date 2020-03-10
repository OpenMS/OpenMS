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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DGoToDialog.h>
#include <ui_Spectrum1DGoToDialog.h>

#include <QtWidgets/QLineEdit>

using namespace std;

namespace OpenMS
{

  Spectrum1DGoToDialog::Spectrum1DGoToDialog(QWidget * parent) :
    QDialog(parent),
    ui_(new Ui::Spectrum1DGoToDialogTemplate)
  {
    ui_->setupUi(this);
  }

  Spectrum1DGoToDialog::~Spectrum1DGoToDialog()
  {
    delete ui_;
  }


  void Spectrum1DGoToDialog::setRange(float min, float max)
  {
    ui_->min_->setText(QString::number(min));
    ui_->max_->setText(QString::number(max));
  }

  void Spectrum1DGoToDialog::setMinMaxOfRange(float min, float max)
  {
    ui_->min_const_->setText(QString("min: ") + QString::number(min));
    ui_->max_const_->setText(QString("max: ") + QString::number(max));
  }

  bool Spectrum1DGoToDialog::checked()
  {
    return ui_->clip_checkbox->checkState() == Qt::Checked;
  }

  void Spectrum1DGoToDialog::fixRange()
  {
    // load from GUI
    float min_mz = ui_->min_->text().toFloat();
    float max_mz = ui_->max_->text().toFloat();

    // ensure correct order of min and max
    if (min_mz > max_mz) swap(min_mz, max_mz);

    // do not allow range of 0 --> extend to 1
    if (min_mz == max_mz)
    {
      min_mz -= 0.5;
      max_mz += 0.5;
    }

    // store in GUI
    ui_->min_->setText(QString::number(min_mz));
    ui_->max_->setText(QString::number(max_mz));
  }

  float Spectrum1DGoToDialog::getMin() const
  {
    return ui_->min_->text().toFloat();
  }

  float Spectrum1DGoToDialog::getMax() const
  {
    return ui_->max_->text().toFloat();
  }

} //namespace OpenMS

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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/PythonSelector.h>
#include <ui_PythonSelector.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/PythonInfo.h>

#include <QString>
#include <QtWidgets/QFileDialog>
#include <QMessageBox>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    PythonSelector::PythonSelector(QWidget* parent) :
    QWidget(parent),
    ui_(new Ui::PythonSelector)
    {
      ui_->setupUi(this);
      
      connect(ui_->btn_browse, SIGNAL(clicked()), this, SLOT(showFileDialog_()));
      connect(ui_->line_edit, SIGNAL(editingFinished()), this, SLOT(validate_()));

      // load/update UI
      ui_->line_edit->setText(last_known_python_exe_.toQString());

      // internally check
      validate_();
    }

    PythonSelector::~PythonSelector()
    {
      delete ui_;
      // TODO: store UI to INI?
    }

    
    void PythonSelector::showFileDialog_()
    {
      QString file_name = QFileDialog::getOpenFileName(this, tr("Specify Python executable"), tr(""), tr(/*valid formats*/ ""));
      if (!file_name.isEmpty())
      {
        ui_->line_edit->setText(file_name); // will not trigger the validator
        emit ui_->line_edit->editingFinished(); // simulate loosing focus or pressing return (to trigger validate_())
      }
    }

    void PythonSelector::validate_()
    {
      String exe = ui_->line_edit->text();
      
      String error;
      bool success = PythonInfo::canRun(exe, error);
      if (success)
      {
        last_known_python_exe_ = exe;
        ui_->label->setText(PythonInfo::getVersion(exe).toQString());
        currently_valid_ = true;
      }
      else
      {
        QMessageBox::warning(0, QString("Python not found"), error.toQString());
        // no need to currently_valid_=false, since we will revert to 'last_known_python_exe_'
      }

      // reset to last known
      ui_->line_edit->setText(last_known_python_exe_.toQString());

      emit valueChanged(last_known_python_exe_.toQString(), currently_valid_);
    }


  }   //namespace Internal
} //namspace OpenMS

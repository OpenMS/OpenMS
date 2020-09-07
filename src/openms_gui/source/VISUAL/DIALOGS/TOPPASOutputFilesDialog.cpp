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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFilesDialog.h>
#include <ui_TOPPASOutputFilesDialog.h>

#include <OpenMS/SYSTEM/File.h>

#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QCompleter>
#include <QtWidgets/QDirModel>
#include <QtCore/QDir>
#include <QtCore/QFileInfo>


namespace OpenMS
{
  TOPPASOutputFilesDialog::TOPPASOutputFilesDialog(const QString& dir_name, int num_jobs)
    : ui_(new Ui::TOPPASOutputFilesDialogTemplate)
  {
    ui_->setupUi(this);
    if (dir_name != "")
    {
      ui_->out_dir->setDirectory(dir_name);
    }
    else
    {
      ui_->out_dir->setDirectory(QDir::currentPath());
    }
    if (num_jobs >= 1)
    {
      ui_->num_jobs_box->setValue(num_jobs);
    }
    
    connect(ui_->ok_button, SIGNAL(clicked()), this, SLOT(checkValidity_()));
    connect(ui_->cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
    
    // make Ok the default (just pressing Enter will run the workflow)
    ui_->ok_button->setFocus();
  }

  TOPPASOutputFilesDialog::~TOPPASOutputFilesDialog()
  {
    delete ui_;
  }

  void TOPPASOutputFilesDialog::showFileDialog()
  {
    ui_->out_dir->showFileDialog();
  }

  QString TOPPASOutputFilesDialog::getDirectory() const
  {
    return ui_->out_dir->getDirectory();
  }

  int TOPPASOutputFilesDialog::getNumJobs() const
  {
    return ui_->num_jobs_box->value();
  }

  void TOPPASOutputFilesDialog::checkValidity_()
  {
    if (!ui_->out_dir->dirNameValid())
    {
      QMessageBox::warning(nullptr, "Invalid directory", "Either the specified path is no directory, or you have no permission to write there.");
      return;
    }

    accept();
  }


} // namespace

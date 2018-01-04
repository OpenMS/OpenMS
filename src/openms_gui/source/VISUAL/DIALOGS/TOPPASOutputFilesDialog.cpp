// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/SYSTEM/File.h>

#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>
#include <QtGui/QCompleter>
#include <QtGui/QDirModel>
#include <QtCore/QDir>
#include <QtCore/QFileInfo>

#include <iostream>

namespace OpenMS
{
  TOPPASOutputFilesDialog::TOPPASOutputFilesDialog(const QString & dir_name, int num_jobs)
  {
    setupUi(this);
    if (dir_name != "")
    {
      line_edit->setText(dir_name);
    }
    else
    {
      line_edit->setText(QDir::currentPath());
    }
    if (num_jobs >= 1)
    {
      num_jobs_box->setValue(num_jobs);
    }
    QCompleter * completer = new QCompleter(this);
    QDirModel * dir_model = new QDirModel(completer);
    dir_model->setFilter(QDir::AllDirs);
    completer->setModel(dir_model);
    line_edit->setCompleter(completer);
    connect(browse_button, SIGNAL(clicked()), this, SLOT(showFileDialog()));
    connect(ok_button, SIGNAL(clicked()), this, SLOT(checkValidity_()));
    connect(cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
    
    // make Ok the default (just pressing Enter will run the workflow)
    ok_button->setFocus();
  }

  void TOPPASOutputFilesDialog::showFileDialog()
  {
    QString dir = File::exists(File::path(line_edit->text())) ? File::path(line_edit->text()).toQString() : "";
    QString selected_dir = QFileDialog::getExistingDirectory(this, tr("Select output directory"), dir);
    if (selected_dir != "")
    {
      line_edit->setText(selected_dir);
    }
  }

  QString TOPPASOutputFilesDialog::getDirectory()
  {
    return line_edit->text();
  }

  int TOPPASOutputFilesDialog::getNumJobs()
  {
    return num_jobs_box->value();
  }

  void TOPPASOutputFilesDialog::checkValidity_()
  {
    if (!dirNameValid(line_edit->text()))
    {
      QMessageBox::warning(nullptr, "Invalid directory", "Either the specified path is no directory, or you have no permission to write there.");
      return;
    }

    accept();
  }

  bool TOPPASOutputFilesDialog::dirNameValid(const QString & dir_name)
  {
    QFileInfo fi(dir_name);
    QString file_name = dir_name;
    if (!file_name.endsWith(QDir::separator()))
    {
      file_name += QDir::separator();
    }
    file_name += "test_file";
    return fi.isDir() && File::writable(file_name);
  }

} // namespace

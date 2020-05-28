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

// OpenMS includes
#include <OpenMS/VISUAL/OutputDirectory.h>
#include <ui_OutputDirectory.h>


#include <OpenMS/SYSTEM/File.h>

#include <QtWidgets/QMessageBox>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QCompleter>
#include <QtWidgets/QDirModel>


namespace OpenMS
{
  OutputDirectory::OutputDirectory(QWidget* parent)
    : QWidget(parent),
      ui_(new Ui::OutputDirectoryTemplate)
  {
    ui_->setupUi(this);
    QCompleter* completer = new QCompleter(this);
    QDirModel* dir_model = new QDirModel(completer);
    dir_model->setFilter(QDir::AllDirs);
    completer->setModel(dir_model);
    ui_->line_edit->setCompleter(completer);

    connect(ui_->browse_button, &QPushButton::clicked, this, &OutputDirectory::showFileDialog);
    connect(ui_->line_edit, &QLineEdit::textChanged, this, &OutputDirectory::textEditChanged_);
  }

  OutputDirectory::~OutputDirectory()
  {
    delete ui_;
  }

  void OutputDirectory::setDirectory(const QString& dir)
  {
    ui_->line_edit->setText(dir);
    emit directoryChanged(dir);
  }

  QString OutputDirectory::getDirectory() const
  {
    return ui_->line_edit->text();
  }

  void OutputDirectory::showFileDialog()
  {
    QString dir = File::exists(File::path(getDirectory())) ? File::path(getDirectory()).toQString() : "";
    QString selected_dir = QFileDialog::getExistingDirectory(this, tr("Select output directory"), dir);
    if (!selected_dir.isEmpty())
    {
      setDirectory(selected_dir); // emits directoryChanged()
    }
  }
  
  void OutputDirectory::textEditChanged_(const QString& /*new_text*/)
  {
    emit directoryChanged(getDirectory());
  }
  

  bool OutputDirectory::dirNameValid() const
  {
    if (!QFileInfo(getDirectory()).isDir()) return false;

    QString file_name = getDirectory();
    if (!file_name.endsWith(QDir::separator()))
    {
      file_name += QDir::separator();
    }
    file_name += "test_file";
    return File::writable(file_name);
  }


} // namespace

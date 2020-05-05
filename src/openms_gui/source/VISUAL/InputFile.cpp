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
#include <OpenMS/VISUAL/InputFile.h>
#include <ui_InputFile.h>

#include <QtWidgets/QMessageBox>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QCompleter>
#include <QtWidgets/QDirModel>


namespace OpenMS
{
  InputFile::InputFile(const QWidget* parent)
    : file_format_filter_(),
      ui_(new Ui::InputFileTemplate)
  {
    ui_->setupUi(this);
    // disable completer for windows, causes crashes
//#ifndef OPENMS_WINDOWSPLATFORM
    QCompleter* completer = new QCompleter(this);
    completer->setModel(new QDirModel(completer));
    ui_->line_edit->setCompleter(completer);
//#endif
    connect(ui_->browse_button, SIGNAL(clicked()), this, SLOT(showFileDialog()));
  }

  InputFile::~InputFile()
  {
    delete ui_;
  }

  void InputFile::setFilename(const QString& filename)
  {
    ui_->line_edit->setText(filename);
  }

  QString InputFile::getFilename() const
  {
    return ui_->line_edit->text();
  }

  void InputFile::setFileFormatFilter(const QString& fff)
  {
    file_format_filter_ = fff;
  }

  void InputFile::showFileDialog()
  {
    QFileInfo fi(getFilename()); // get path from current file as starting directory for selection
    
    QString file_name = QFileDialog::getOpenFileName(this, tr("Specify input file"), fi.path(), file_format_filter_);
    if (!file_name.isEmpty())
    {
      ui_->line_edit->setText(file_name);
    }
  }


} // namespace

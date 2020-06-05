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
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>
#include <ui_TOPPASInputFilesDialog.h>

#include <OpenMS/VISUAL/InputFileList.h>

namespace OpenMS
{
  TOPPASInputFilesDialog::TOPPASInputFilesDialog(const QStringList& list, const QString& cwd, QWidget* parent)
    : QDialog(parent),
      ui_(new Ui::TOPPASInputFilesDialogTemplate)
  {
    ui_->setupUi(this);
    ifl_ = (InputFileList*)ui_->input_file_list;
    ifl_->setCWD(cwd);
    ifl_->setFilenames(list);

    connect(ui_->ok_button, SIGNAL(clicked()), this, SLOT(accept()));
    connect(ui_->cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
    setAcceptDrops(true);
  }

  TOPPASInputFilesDialog::~TOPPASInputFilesDialog()
  {
    delete ui_;
  }

  void TOPPASInputFilesDialog::getFilenames(QStringList& files) const
  {
    ifl_->getFilenames(files);
    if (ui_->flag_sort_list->isChecked())
      files.sort();
  }

  const QString& TOPPASInputFilesDialog::getCWD() const
  {
    return ifl_->getCWD();
  }



} // namespace

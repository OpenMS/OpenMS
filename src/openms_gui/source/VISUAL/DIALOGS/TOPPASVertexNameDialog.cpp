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
#include <OpenMS/VISUAL/DIALOGS/TOPPASVertexNameDialog.h>
#include <ui_TOPPASVertexNameDialog.h>

#include <QRegExpValidator> 

#include <iostream>

namespace OpenMS
{
  TOPPASVertexNameDialog::TOPPASVertexNameDialog(const QString& name, const QString& input_regex)
    : ui_(new Ui::TOPPASVertexNameDialogTemplate)
  {
    ui_->setupUi(this);
    
    if (!input_regex.isEmpty())
    {
      QRegExp rx(input_regex);
      QRegExpValidator* v = new QRegExpValidator(rx, ui_->line_edit);
      ui_->line_edit->setValidator(v);
    }

    ui_->line_edit->setText(name);
    connect(ui_->ok_button, SIGNAL(clicked()), this, SLOT(accept()));
    connect(ui_->cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
  }

  TOPPASVertexNameDialog::~TOPPASVertexNameDialog()
  {
    delete ui_;
  }

  QString TOPPASVertexNameDialog::getName()
  {
    
    return ui_->line_edit->text();
  }

} // namespace

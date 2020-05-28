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

#include <OpenMS/VISUAL/DIALOGS/PythonModuleRequirement.h>
#include <ui_PythonModuleRequirement.h>

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
    PythonModuleRequirement::PythonModuleRequirement(QWidget* parent) :
      QWidget(parent),
      ui_(new Ui::PythonModuleRequirement)
    {
      ui_->setupUi(this);
    }

    // slot
    void PythonModuleRequirement::validate(const QString& python_exe)
    {
      QStringList valid_modules;
      QStringList missing_modules;
      ui_->lbl_modules->setText(" ... updating ... ");
      for (const auto& s : required_modules_)
      {
        if (PythonInfo::isPackageInstalled(python_exe, s)) valid_modules.push_back(s);
        else missing_modules.push_back(s);
      }
      emit valueChanged(valid_modules, missing_modules);
      QString text = "<ul>";
      if (!valid_modules.empty()) text += QString("<li> [<code style = \"color: green\">%1</code>] present").arg(valid_modules.join(", "));
      if (!missing_modules.empty()) text += QString("<li> [<code style = \"color: red\">%1</code>] missing").arg(missing_modules.join(", "));
      text += "</ul>";
      ui_->lbl_modules->setText(text);
      // if no modules are missing, we are good to go...
      is_ready_ = missing_modules.empty();
    }

    PythonModuleRequirement::~PythonModuleRequirement()
    {
      delete ui_;
      // TODO: store UI to INI?
    }

    void PythonModuleRequirement::setTitle(const QString& title)
    {
      ui_->box->setTitle(title);
    }

    void PythonModuleRequirement::setRequiredModules(const QStringList& m)
    {
      required_modules_ = m;
    }

    void PythonModuleRequirement::setFreeText(const QString& text)
    {
      ui_->lbl_freetext->setText(text);
    }


  } //namespace Internal
} //namspace OpenMS


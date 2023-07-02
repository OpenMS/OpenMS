// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/VISUAL/APPLICATIONS/FLASHDeconvWizardBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/VISUAL/DIALOGS/FLASHDeconvTabWidget.h>
#include <ui_FLASHDeconvWizardBase.h>

// Qt
#include <QDesktopServices>
#include <QMessageBox>
#include <QSettings>

using namespace std;
using namespace OpenMS;

namespace OpenMS
{
  using namespace Internal;

  FLASHDeconvWizardBase::FLASHDeconvWizardBase(QWidget* parent) :
      QMainWindow(parent), DefaultParamHandler("FLASHDeconvWizardBase"),
      ui(new Ui::FLASHDeconvWizardBase)
  {
    ui->setupUi(this);
    QSettings settings("OpenMS", "FLASHDeconvWizard");
    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());
    setWindowTitle("FLASHDeconvWizard");
    setWindowIcon(QIcon(":/FLASHDeconvWizard.png"));

    FLASHDeconvTabWidget* cwidget = new FLASHDeconvTabWidget(this);
    setCentralWidget(cwidget);
  }


  FLASHDeconvWizardBase::~FLASHDeconvWizardBase()
  {
    delete ui;
  }


  void FLASHDeconvWizardBase::showAboutDialog()
  {
    QApplicationTOPP::showAboutDialog(this, "FLASHDeconvWizard");
  }

  void OpenMS::FLASHDeconvWizardBase::on_actionExit_triggered()
  {
    QApplicationTOPP::exit();
  }

  void OpenMS::FLASHDeconvWizardBase::on_actionVisit_FLASHDeconv_homepage_triggered()
  {
    const char* url = "https://openms.de/application/flashdeconv/";
    if (!QDesktopServices::openUrl(QUrl(url)))
    {
      QMessageBox::warning(nullptr, "Cannot open browser. Please check your default browser settings.", QString(url));
    }
  }

  void OpenMS::FLASHDeconvWizardBase::on_actionReport_new_issue_triggered()
  {
    const char* url = "https://github.com/OpenMS/OpenMS/issues";
    if (!QDesktopServices::openUrl(QUrl(url)))
    {
      QMessageBox::warning(nullptr, "Cannot open browser. Please check your default browser settings.", QString(url));
    }
  }

} // namespace OpenMS

// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/VISUAL/APPLICATIONS/FLASHQuantWizardBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/VISUAL/DIALOGS/FLASHQuantTabWidget.h>
#include <ui_FLASHQuantWizardBase.h>

// Qt
#include <QDesktopServices>
#include <QMessageBox>
#include <QSettings>

using namespace std;
using namespace OpenMS;

namespace OpenMS
{
  using namespace Internal;

  FLASHQuantWizardBase::FLASHQuantWizardBase(QWidget* parent) :
    QMainWindow(parent),
    DefaultParamHandler("FLASHQuantWizardBase"),
    //clipboard_scene_(nullptr),
    ui(new Ui::FLASHQuantWizardBase)
  {
    ui->setupUi(this);
    QSettings settings("OpenMS", "FLASHQuantWizard");
    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());
    setWindowTitle("FLASHQuantWizard");
    setWindowIcon(QIcon(":/FLASHDeconvWizard.png"));

    FLASHQuantTabWidget* cwidget = new FLASHQuantTabWidget(this);
    setCentralWidget(cwidget);
  }


  FLASHQuantWizardBase::~FLASHQuantWizardBase()
  {
    delete ui;
  }


  void FLASHQuantWizardBase::showAboutDialog()
  {
    QApplicationTOPP::showAboutDialog(this, "FLASHQuantWizard");
  }

  void OpenMS::FLASHQuantWizardBase::on_actionExit_triggered()
  {
      QApplicationTOPP::exit();
  }

  void OpenMS::FLASHQuantWizardBase::on_actionVisit_FLASHQuant_homepage_triggered()
  {
    const char* url = "https://www.openms.de/comp/FLASHQuant/";
    if (!QDesktopServices::openUrl(QUrl(url)))
    {
      QMessageBox::warning(0, "Cannot open browser. Please check your default browser settings.", QString(url));
    }
  }

  void OpenMS::FLASHQuantWizardBase::on_actionReport_new_issue_triggered()
  {
    const char* url = "https://github.com/OpenMS/OpenMS/issues"; // TODO: change?
    if (!QDesktopServices::openUrl(QUrl(url)))
    {
      QMessageBox::warning(0, "Cannot open browser. Please check your default browser settings.", QString(url));
    }
  }

} //namespace OpenMS






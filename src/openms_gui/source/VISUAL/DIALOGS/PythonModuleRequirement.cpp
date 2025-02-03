// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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


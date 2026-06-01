// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class PythonModuleRequirement;
}

namespace OpenMS
{
  namespace Internal
  {
    /// Given a list of python modules which are required, this widget checks them and
    /// displays the current status
    class OPENMS_GUI_DLLAPI PythonModuleRequirement : public QWidget
    {
      Q_OBJECT
      
    public:
      explicit PythonModuleRequirement(QWidget* parent = nullptr);
      ~PythonModuleRequirement();

      /// change the label of the surrounding box
      void setTitle(const QString& title);

      /// a list of python modules required for a certain functionality/script
      void setRequiredModules(const QStringList& m);

      /// some arbitrary description for the user to display statically
      void setFreeText(const QString& text);

      /// are all modules present?
      bool isReady() const { return is_ready_;};


    signals:
      /// emitted whenever the requirement check was executed...
      void valueChanged(QStringList& valid_modules, QStringList& missing_modules);


    public slots:
      /// re-evaluate the presence of modules, based on a new python version
      void validate(const QString& python_exe);

    private:
      QStringList required_modules_; ///< list of modules which are needed (order might be important -- know your Python...)
      QString info_text_; ///< additional text to display for the user
      bool is_ready_ = false; ///< all modules are present and the app is good to go

      Ui::PythonModuleRequirement* ui_;
    };

  } // ns Internal
} // ns OpenMS

// this is required to allow Ui_SwathTabWidget (auto UIC'd from .ui) to have a PythonModuleRequirement member
using PythonModuleRequirement = OpenMS::Internal::PythonModuleRequirement;

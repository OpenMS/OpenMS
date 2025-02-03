// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class PythonSelector;
}

namespace OpenMS
{
  namespace Internal
  {
    /// A QLineEdit + Browse button to have the user select a local python installation
    /// By default, 'python' is used
    class OPENMS_GUI_DLLAPI PythonSelector : public QWidget
    {
      Q_OBJECT

    public:
      explicit PythonSelector(QWidget* parent = nullptr);
      ~PythonSelector();

      const String& getLastPython() const
      {
        return last_known_python_exe_;
      }

    signals:
      /// emitted whenever the line-edit has new values for the current python executable
      /// @param last_known_python_exe The currently best guess where python can be found
      /// @param valid_python Is the python executable given in @p last_known_python_exe callable?
      void valueChanged(QString last_known_python_exe, bool valid_python);
      
    
    private slots:
      void showFileDialog_();

      void validate_();

    private:
      String last_known_python_exe_ = "python"; ///< initial guess or last valid user input
      bool currently_valid_ = false; ///< unless proven otherwise by 'validate_()'

      Ui::PythonSelector* ui_;
    };

  }
} // ns OpenMS

// this is required to allow Ui_SwathTabWidget (auto UIC'd from .ui) to have a PythonSelector member
using PythonSelector = OpenMS::Internal::PythonSelector;

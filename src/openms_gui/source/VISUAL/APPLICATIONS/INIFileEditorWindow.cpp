// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/APPLICATIONS/INIFileEditorWindow.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtWidgets/QToolBar>
#include <QtCore/QString>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QCheckBox>
#include <QCloseEvent>

using namespace std;

namespace OpenMS
{

  INIFileEditorWindow::INIFileEditorWindow(QWidget* parent) :
    QMainWindow(parent),
    current_path_(".")
  {
    setWindowTitle("INIFileEditor");
    setWindowIcon(QIcon(":/INIFileEditor.png"));

    //create central widget and layout
    QWidget* central_widget = new QWidget;
    setCentralWidget(central_widget);
    QGridLayout* layout = new QGridLayout(central_widget);

    //create advanced check box and ParamEditor and connect them
    editor_ = new ParamEditor(central_widget);
    layout->addWidget(editor_, 0, 0, 1, 2);

    QMenu* file = new QMenu("&File", this);
    menuBar()->addMenu(file);
    file->addAction("&Open", this, SLOT(openFile()), Qt::CTRL + Qt::Key_O);
    file->addSeparator();
    file->addAction("&Save", this, SLOT(saveFile()), Qt::CTRL + Qt::Key_S);
    file->addAction("Save &As", this, SLOT(saveFileAs()));
    file->addSeparator();
    file->addAction("&Quit", this, SLOT(close()));

    // we connect the "changes state"(changes made/no changes) signal from the ParamEditor to the window title updating slot
    connect(editor_, SIGNAL(modified(bool)), this, SLOT(updateWindowTitle(bool)));

    setMinimumSize(600, 600);
  }

  bool INIFileEditorWindow::openFile(const String& filename)
  {
    if (filename.empty())
    {
      filename_ = QFileDialog::getOpenFileName(this, tr("Open ini file"), current_path_.toQString(), tr("ini files (*.ini);; all files (*.*)"));
    }
    else
    {
      filename_ = filename.c_str();
    }

    if (!filename_.isEmpty())
    {
      if (File::readable(filename_.toStdString()))
      {

        param_.clear();
        ParamXMLFile paramFile;
        try
        {
          paramFile.load(filename_.toStdString(), param_);
          editor_->load(param_);
          updateWindowTitle(editor_->isModified());
          return true;
        }
        catch (Exception::BaseException& e)
        {
          OPENMS_LOG_ERROR << "Error while parsing file '" << filename_.toStdString() << "'\n";
          OPENMS_LOG_ERROR << e << "\n";
        }
      }

      QMessageBox::critical(this, "Error opening file", ("The file '" + filename_.toStdString() + "' does not exist, is not readable or not a proper INI file!").c_str());
    }
    return false;
  }

  bool INIFileEditorWindow::saveFile()
  {
    if (filename_.isEmpty())
    {
      return false;
    }

    editor_->store();

    ParamXMLFile paramFile;
    paramFile.store(filename_.toStdString(), param_);
    updateWindowTitle(editor_->isModified());
    return true;
  }

  bool INIFileEditorWindow::saveFileAs()
  {
    filename_ = QFileDialog::getSaveFileName(this, tr("Save ini file"), current_path_.toQString(), tr("ini files (*.ini)"));
    if (!filename_.isEmpty())
    {
      if (!filename_.endsWith(".ini"))
        filename_.append(".ini");

      editor_->store();

      ParamXMLFile paramFile;
      paramFile.store(filename_.toStdString(), param_);
      updateWindowTitle(editor_->isModified());
      return true;
    }
    return false;
  }

  void INIFileEditorWindow::closeEvent(QCloseEvent* event)
  {
    if (editor_->isModified())
    {
      QMessageBox::StandardButton result = QMessageBox::question(this, "Save?", "Do you want to save your changes?", QMessageBox::Ok | QMessageBox::Cancel | QMessageBox::Discard);
      if (result == QMessageBox::Ok)
      {
        if (saveFile())
        {
          event->accept();
        }
        else
        {
          event->ignore();
        }
      }
      else if (result == QMessageBox::Cancel)
      {
        event->ignore();
      }
      else
      {
        event->accept();
      }
    }
    else
    {
      event->accept();
    }
  }

  void INIFileEditorWindow::updateWindowTitle(bool update)
  {
    //update window title
    if (update)
    {
      setWindowTitle((File::basename(filename_) + " * - INIFileEditor").toQString());
    }
    else
    {
      setWindowTitle((File::basename(filename_) + " - INIFileEditor").toQString());
    }

    //update last path as well
    current_path_ = File::path(filename_);
  }

}

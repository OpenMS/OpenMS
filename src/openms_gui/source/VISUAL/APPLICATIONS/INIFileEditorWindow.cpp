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
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

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
    file->addAction("&Open", this, [this]() { openFile(); })->setShortcut(Qt::CTRL | Qt::Key_O);
    file->addSeparator();
    file->addAction("&Save", this, &INIFileEditorWindow::saveFile)->setShortcut(Qt::CTRL | Qt::Key_S);
    file->addAction("Save &As", this, &INIFileEditorWindow::saveFileAs);
    file->addSeparator();
    file->addAction("&Quit", this, &QWidget::close);

    // we connect the "changes state"(changes made/no changes) signal from the ParamEditor to the window title updating slot
    connect(editor_, &ParamEditor::modified, this, &INIFileEditorWindow::updateWindowTitle);

    setMinimumSize(600, 600);
  }

  bool INIFileEditorWindow::openFile(const std::string& filename)
  {
    if (filename.empty())
    {
      filename_ = QFileDialog::getOpenFileName(this, tr("Open ini file"), toQString(current_path_), tr("ini files (*.ini);; all files (*.*)"));
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

    // refuse to write rather than persist the outdated param when an edit could not be committed
    // (e.g. the open editor holds a value that violates the parameter's restrictions)
    if (!editor_->store())
    {
      QMessageBox::warning(this, "Not saved", "The file was not saved because the value currently being edited could not be applied. Please correct it and save again.");
      return false;
    }

    ParamXMLFile paramFile;
    try
    {
      paramFile.store(filename_.toStdString(), param_);
    }
    catch (Exception::BaseException& e)
    {
      // store() now refuses non-atomic writes (read-only target, directory, rename failure) by
      // throwing. The document is left as-is: editor_->store() commits the edit into the model but
      // no longer marks it clean, so any unsaved changes are still reflected and closing will still
      // prompt to save. If there was nothing to save, it correctly stays clean.
      QMessageBox::critical(this, "Error saving file", ("The file '" + filename_.toStdString() + "' could not be saved: " + e.getMessage()).c_str());
      return false;
    }
    editor_->setModified(false); // persisted successfully -- the document is now clean
    updateWindowTitle(editor_->isModified());
    return true;
  }

  bool INIFileEditorWindow::saveFileAs()
  {
    QString filename = QFileDialog::getSaveFileName(this, tr("Save ini file"), toQString(current_path_), tr("ini files (*.ini)"));
    if (filename.isEmpty())
    {
      return false;
    }

    if (!filename.endsWith(".ini"))
      filename.append(".ini");

    if (!editor_->store())
    {
      QMessageBox::warning(this, "Not saved", "The file was not saved because the value currently being edited could not be applied. Please correct it and save again.");
      return false;
    }

    ParamXMLFile paramFile;
    try
    {
      paramFile.store(filename.toStdString(), param_);
    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::critical(this, "Error saving file", ("The file '" + filename.toStdString() + "' could not be saved: " + e.getMessage()).c_str());
      return false;
    }
    // adopt the new name only once the write has actually succeeded, so a cancelled or failed
    // Save As does not replace the name of the file currently being edited
    filename_ = filename;
    editor_->setModified(false); // persisted successfully -- the document is now clean
    updateWindowTitle(editor_->isModified());
    return true;
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
      setWindowTitle(toQString((File::basename(fromQString(filename_)) + " * - INIFileEditor")));
    }
    else
    {
      setWindowTitle(toQString((File::basename(fromQString(filename_)) + " - INIFileEditor")));
    }

    //update last path as well
    current_path_ = File::path(fromQString(filename_));
  }

}

// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
    if (filename == "")
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

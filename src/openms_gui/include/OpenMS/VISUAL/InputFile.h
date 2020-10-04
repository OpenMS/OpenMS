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

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class InputFileTemplate;
}

namespace OpenMS
{
  /**
      @brief A simple widget with a line-edit and a browse button to choose filenames

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI InputFile :
    public QWidget
  {
    Q_OBJECT
  public:
    /// Constructor
    InputFile(QWidget* parent);
    /// Destructor
    ~InputFile();

    /// support drag'n'drop of files from OS window manager
    void dragEnterEvent(QDragEnterEvent* e) override;
    /// support drag'n'drop of files from OS window manager
    void dropEvent(QDropEvent* e) override;
    void dragMoveEvent(QDragMoveEvent* pEvent) override;

    /// Sets the text in the line-edit
    void setFilename(const QString& filename);

    /// Returns the filename currently set in the line-edit
    QString getFilename() const;

    /// Users can only choose certain filetypes, e.g. "Transition sqLite file (*.pqp)"
    void setFileFormatFilter(const QString& fff);

    /// get the CWD (according to most recently added file)
    const QString& getCWD() const;
    /// set the current working directory (for opening files). If the input is not empty, the cwd will not be altered, unless @p force is used
    void setCWD(const QString& cwd, bool force = false);
 
  signals:
    /// emitted when a new file is added (by drag'n'drop or 'Browse' button)
    void updatedCWD(QString new_cwd);
    /// emitted when a new file is added (by drag'n'drop or 'Browse' button)
    void updatedFile(QString new_path);

  public slots:

    /// Lets the user select the file via a file dialog
    void showFileDialog();


  protected:
    /// optional filter during file browsing
    QString file_format_filter_;
    /// the current working directory according to the last file added
    QString cwd_;

  private:
    Ui::InputFileTemplate* ui_;
  };

}

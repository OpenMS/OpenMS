// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_VISUAL_DIALOGS_SAVEIMAGEDIALOG_H
#define OPENMS_VISUAL_DIALOGS_SAVEIMAGEDIALOG_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtGui/QDialog>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>
#include <QtGui/QCheckBox>

namespace OpenMS
{
  /**
      @brief Dialog for saving an image.

      @image html SaveImageDialog.png

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI SaveImageDialog :
    public QDialog
  {
    Q_OBJECT

public:
    ///Constructor
    SaveImageDialog(QWidget * parent = nullptr);
    ///set size and size ratio
    void setSize(int x, int y);
    ///accessors for the width
    int getXSize();
    ///accessors for the height
    int getYSize();
    ///accessors for the format
    QString getFormat();

public slots:
    ///changes width keeping proportions
    void xSizeChanged(const QString & s);
    ///changes height keeping proportions
    void ySizeChanged(const QString & s);
    ///set size ratio when proportions checkbox is activated
    void proportionsActivated(bool state);
    ///checks if the values for width and height are ok before accepting the dialog
    void checkSize();

private:
    //format
    QComboBox * format_;
    //size
    QLineEdit * size_x_;
    QLineEdit * size_y_;
    QCheckBox * size_proportions_;
    //ratio size_x_/size_y_
    float size_ratio_;

    //set the size ratio (width/height)
    void setSizeRatio_(float r);
  };
}
#endif

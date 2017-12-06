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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_MULTIGRADIENTSELECTOR_H
#define OPENMS_VISUAL_MULTIGRADIENTSELECTOR_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/VISUAL/MultiGradient.h>

//QT
#include <QtGui/QWidget>
class QPaintEvent;
class QMouseEvent;
class QKeyEvent;
class QContextMenuEvent;

//std lib
#include <vector>

namespace OpenMS
{

  /**
      @brief A widget witch allows constructing gradients of multiple colors.

      @image html MultiGradientSelector.png

      The above example image shows a MultiGradientSelector.

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI MultiGradientSelector :
    public QWidget
  {
    Q_OBJECT
public:
    ///Constructor
    MultiGradientSelector(QWidget * parent = nullptr);
    ///Destructor
    ~MultiGradientSelector() override;

    ///returns a const reference to the gradient
    const MultiGradient & gradient() const;
    ///returns a mutable reference to the gradient
    MultiGradient & gradient();

    /// sets the interpolation mode
    void setInterpolationMode(MultiGradient::InterpolationMode mode);
    /// returns the interpolation mode
    MultiGradient::InterpolationMode getInterpolationMode() const;

public slots:
    /// sets what interpolation mode is used
    void stairsInterpolation(bool state);

protected:
    ///@name re-implemented Qt events
    //@{
    void paintEvent(QPaintEvent * e) override;
    void mousePressEvent(QMouseEvent * e) override;
    void mouseMoveEvent(QMouseEvent * e) override;
    void mouseReleaseEvent(QMouseEvent * e) override;
    void mouseDoubleClickEvent(QMouseEvent * e) override;
    void keyPressEvent(QKeyEvent * e) override;
    void contextMenuEvent(QContextMenuEvent * e) override;
    //@}

    // the actual gradient
    MultiGradient gradient_;

    // border margin
    Int margin_;
    // height of the gradient area
    Int gradient_area_width_;
    // height of the lever area
    Int lever_area_height_;

    //position (0-100) in the vector of the selected lever
    Int selected_;
    //color of the selected lever
    QColor selected_color_;

    //stores if the left button is pressed
    bool left_button_pressed_;
  };

}
#endif // OPENMS_VISUAL_MULTIGRADIENTSELECTOR_H

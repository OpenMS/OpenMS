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
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_VISUALIZER_METAINFOVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_METAINFOVISUALIZER_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

//STL
#include <vector>
#include <utility>

class QAbstractButton;
class QButtonGroup;

namespace OpenMS
{
  /**
      @brief MetaInfoVisualizer is a visualizer class for all classes that use one MetaInfo object as member.

      Meta information is an array of Type-Name-Value tuples. Classes that have a MetaInfo objects as a member can use this class to edit the MetaInfo object.
  */
  class OPENMS_GUI_DLLAPI MetaInfoVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<MetaInfoInterface>
  {
    Q_OBJECT

public:
    ///Constructor
    MetaInfoVisualizer(bool editable = false, QWidget * parent = nullptr);

    //Docu in base class
    void load(MetaInfoInterface & m);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    /// Adds a new Type-Value pair to the MetaInfo Object.
    void add_();
    /// Clears out all fields.
    void clear_();
    /// Removes a selected Type-Value pair from the MetaInfo Object.
    void remove_(int);
    ///Undo the changes made in the GUI.
    void undo_();

protected:
    /// Loads all Type-Value pairs one after another.
    void loadData_(UInt index);

    /** @name Edit fields for new Type-Value pair.
      */
    //@{
    QLineEdit * newkey_;
    QLineEdit * newvalue_;
    QLineEdit * newdescription_;
    //@}

    ///@name Arrays of pointers to objects for temporary metaInfo data
    //@{
    std::vector<std::pair<UInt, QLineEdit *> > metainfoptr_;
    std::vector<std::pair<UInt, QLabel *> > metalabels_;
    std::vector<std::pair<UInt, QAbstractButton *> > metabuttons_;
    //@}

    ///@name Edit fields and buttons
    //@{
    QPushButton * addbutton_;
    QPushButton * clearbutton_;
    QButtonGroup * buttongroup_;
    //@}

    /// Counter to keep track of the actual row in the layout.
    int nextrow_;

    /// The layout to display the Type-Value pairs.
    QGridLayout * viewlayout_;

    /// Container for metainfo data.
    std::vector<UInt> keys_;
  };


}
#endif

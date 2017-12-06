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

#ifndef OPENMS_VISUAL_VISUALIZER_GRADIENTVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_GRADIENTVISUALIZER_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>
#include <OpenMS/METADATA/Gradient.h>

//STL
#include <vector>

class QIntValidator;

namespace OpenMS
{
  /**
      @brief GradientVisualizer is a visualizer class for objects of type gradient.

      Each HPLC objects contains a gradient object. A gradient objects contains a list of eluents, timepoints and percentage values. Values can be added to the list, or the whole list can be deleted.
  */
  class OPENMS_GUI_DLLAPI GradientVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<Gradient>
  {
    Q_OBJECT

public:

    ///Constructor
    GradientVisualizer(bool editable = false, QWidget * parent = nullptr);

    //Docu in base class
    void load(Gradient & g);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    /// Add new timepoint to the list
    void addTimepoint_();
    /// Add new eluent to the list
    void addEluent_();
    ///Delete all data from gradient
    void deleteData_();
    ///Undo the changes made in the GUI.
    void undo_();

protected:
    /// Loads a list of eluent, timepoint and percentage triplets.
    void loadData_();
    /// Remove all data from layout
    void removeData_();


    /** @name Edit fields for new eluent-timepoint-percentage-triplets.
      */
    //@{
    QLineEdit * new_eluent_;
    QLineEdit * new_timepoint_;
    //@}

    /** @name Arrays of string values containing eluent, timepoint and percentage values.
    */
    //@{
    std::vector<String> eluents_;
    std::vector<Int> timepoints_;
    //@}

    /** @name Some buttons.
    */
    //@{
    QPushButton * add_eluent_button_;
    QPushButton * add_timepoint_button_;
    QPushButton * removebutton_;
    //@}

    /// Array of temporary pointers to gradient edit fields
    std::vector<QLineEdit *> gradientdata_;

    /// Array of temporary pointers to gradient labels
    std::vector<QLabel *> gradientlabel_;

    /// Pointer to fields with actual data
    QLineEdit * percentage_;

    /// A validator to check the input for the new timepoint.
    QIntValidator * timepoint_vali_;

    /// Counter to keep track of the actual row in the layout.
    int nextrow_;

    /// The layout to display the eluents, timepoints and percentages.
    QGridLayout * viewlayout_;

    //Docu in base class
    void update_() override;
  };


}
#endif

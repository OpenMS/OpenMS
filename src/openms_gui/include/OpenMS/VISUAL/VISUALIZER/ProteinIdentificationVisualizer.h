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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  class MetaDataBrowser;

  /**
      @brief Class that displays all meta information for ProteinIdentification objects

      This class provides all functionality to view the meta information of an object of type Identification.
  */
  class OPENMS_GUI_DLLAPI ProteinIdentificationVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<ProteinIdentification>
  {
    Q_OBJECT

public:
    ///Constructor
    ProteinIdentificationVisualizer(bool editable = false, QWidget * parent = 0, MetaDataBrowser * caller = nullptr);

    /// Loads the meta data from the object to the viewer. Gets the id of the item in the tree as parameter.
    void load(ProteinIdentification & s, int tree_item_id);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    /**
        @brief Updates the tree by calling MetaDataBrowser::updateProteinHits()

        Calls MetaDataBrowser::updateProteinHits().<br>
        Updates the tree depending of the protein significance threshold.<br>
        Only ProteinHits with a score superior or equal to the current threshold will be displayed.
    */
    void updateTree_();

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    /// Pointer to MetaDataBrowser
    MetaDataBrowser * pidv_caller_;
    /// The id of the item in the tree
    int tree_id_;

    ///@name Edit fields and buttons
    //@{
    QLineEdit * engine_;
    QLineEdit * engine_version_;
    QLineEdit * identification_date_;
    QLineEdit * identification_threshold_;
    QLineEdit * identifier_;
    QLineEdit * score_type_;
    QComboBox * higher_better_;

    QLineEdit * db_;
    QLineEdit * db_version_;
    QLineEdit * taxonomy_;
    QLineEdit * charges_;
    QLineEdit * missed_cleavages_;
    QLineEdit * peak_tolerance_;
    QLineEdit * precursor_tolerance_;
    QComboBox * mass_type_;
    QLineEdit * enzyme_;
    //@}

    /// Threshold for filtering by score
    QLineEdit * filter_threshold_;
  };
}

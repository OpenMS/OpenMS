// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ID/MapPeptideIdentifications.h>

#include <functional>
#include <list>

//QT
#include <QtWidgets/QDialog>
#include <QtWidgets/QTreeWidget>
class QTreeWidgetItem;
class QPushButton;
class QStackedWidget;
class QHBoxLayout;
class QVBoxLayout;
class QGridLayout;

namespace OpenMS
{
  class BaseVisualizerGUI;
  class Acquisition;
  class AcquisitionInfo;
  class ContactPerson;
  class ExperimentalSettings;
  class Gradient;
  class HPLC;
  class PeptideIdentification;
  class Instrument;
  class IonDetector;
  class IonSource;
  class MassAnalyzer;
  class MetaInfo;
  class MetaInfoDescription;
  class MetaInfoInterface;
  class MetaInfoRegistry;
  class PeptideHit;
  class Precursor;
  class DataProcessing;
  class ProteinHit;
  class ProteinIdentification;
  class Sample;
  class Software;
  class SourceFile;
  class SpectrumSettings;
  class DocumentIdentifier;
  class Product;

  /**
      @brief A meta data visualization widget

      @image html MetaDataBrowser.png

      It contains a tree view showing all objects of the meta data to be viewed in hierarchical order.

      The meta info data of the tree items are shown in the right part of the viewer, when they are selected in the tree.

      If the data has been modified exec() returns @em true . Otherwise @em false is returned.

      @improvement Add generic mechanism to add items to data vectors e.g. for Instrument - IonSource (Hiwi)

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI MetaDataBrowser :
    public QDialog
  {
    Q_OBJECT

public:

    /// Constructor with flag for edit mode
    MetaDataBrowser(bool editable = false, QWidget * parent = nullptr, bool modal = false);

    /// Adds a peak map
    void add(PeakMap & exp)
    {
      add(static_cast<ExperimentalSettings &>(exp));
      treeview_->expandItem(treeview_->findItems(QString::number(0), Qt::MatchExactly, 1).first());
    }

    /// Adds a peak spectrum
    void add(MSSpectrum & spectrum)
    {
      //spectrum settings
      add(static_cast<SpectrumSettings &>(spectrum));

      //MetaInfoDescriptions
      for (Size i = 0; i < spectrum.getFloatDataArrays().size(); ++i)
      {
        add(spectrum.getFloatDataArrays()[i]);
      }
      for (Size i = 0; i < spectrum.getIntegerDataArrays().size(); ++i)
      {
        add(spectrum.getIntegerDataArrays()[i]);
      }
      for (Size i = 0; i < spectrum.getStringDataArrays().size(); ++i)
      {
        add(spectrum.getStringDataArrays()[i]);
      }

      add(static_cast<MetaInfoInterface &>(spectrum));

      treeview_->expandItem(treeview_->findItems(QString::number(0), Qt::MatchExactly, 1).first());
    }

    /**
      @brief Adds a feature map

      Its identifications are shown (and edited) as peptide and protein identifications (see MapPeptideIdentifications);
      saving writes edited ones back to the map.
    */
    void add(FeatureMap& map);

    /// Adds the feature at @p index of @p map, with its identifications as in add(FeatureMap&)
    void add(FeatureMap& map, Size index);

    /// Adds a consensus map, with its identifications as in add(FeatureMap&)
    void add(ConsensusMap& map);

    /// Adds the consensus feature at @p index of @p map, with its identifications as in add(FeatureMap&)
    void add(ConsensusMap& map, Size index);

    /**
        @brief A generic function to add data.

        The meta data information of all classes that for which a visualize_ method exists can be visualized.
    */
    template <class MetaDataType>
    void add(MetaDataType & meta_data_object)
    {
      visualize_(meta_data_object);
      treeview_->expandItem(treeview_->findItems(QString::number(0), Qt::MatchExactly, 1).first());
    }

    /// Check if mode is editable or not
    bool isEditable() const;

    /// Defines friend classes that can use the functionality of the subclasses.
    friend class ProteinIdentificationVisualizer;
    friend class PeptideIdentificationVisualizer;

public slots:

    /// Set a list of error strings due to invalid date format.
    void setStatus(const std::string& status);

protected slots:

    /// Raises the corresponding viewer from the widget stack according to the item selected in the tree.
    void showDetails_();

    /// Saves all changes and close explorer
    void saveAll_();

protected:

    ///@name Visualizer for the different classes
    //@{
    void visualize_(ExperimentalSettings & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(SpectrumSettings & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(MetaInfoInterface & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(Sample & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(HPLC & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(Gradient & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(Software & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(ScanWindow & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(SourceFile & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(ContactPerson & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(Instrument & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(IonSource & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(IonDetector & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(MassAnalyzer & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(DataProcessingPtr & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(ProteinIdentification & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(ProteinHit & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(PeptideHit & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(Acquisition & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(AcquisitionInfo & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(MetaInfoDescription & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(Precursor & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(Product & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(InstrumentSettings & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(PeptideIdentification & meta, QTreeWidgetItem * parent = nullptr);
    void visualize_(DocumentIdentifier & meta, QTreeWidgetItem * parent = nullptr);
    //@}

    /// Visualizes all elements of a container
    template <typename ContainerType>
    void visualizeAll_(ContainerType & container, QTreeWidgetItem * parent)
    {
      for (typename ContainerType::iterator it = container.begin(); it != container.end(); ++it)
      {
        visualize_(*it, parent);
      }
    }

    /// Connects the Signals of all visualiser classes with Slot setStatus()
    void connectVisualizer_(BaseVisualizerGUI * ptr);

    /// Filters hits according to a score @a threshold. Takes the score orientation into account
    void filterHits_(double threshold, bool higher_better, int tree_item_id);
    /// Shows hits.
    void showAllHits_(int tree_item_id);

    /// The peptide and protein identifications of an added map, to show and edit them (written back by saveAll_())
    template <class MapType>
    MapPeptideIdentifications& addMapIdentifications_(MapType& map);

    /// A list of setting errors due to invalid formats.
    std::string status_list_;

    /// The identifications of the added maps (stable addresses: the visualizers edit them)
    std::list<MapPeptideIdentifications> map_identifications_;
    /// Writes the edited identifications back to their maps
    std::vector<std::function<void()>> write_back_;

    /// Indicates the mode
    bool editable_;

    /// A widgetstack that keeps track of all widgets.
    QStackedWidget * ws_;
    /// Save button
    QPushButton * saveallbutton_;
    /// Close Button
    QPushButton * closebutton_;
    /// Cancel Button
    QPushButton * cancelbutton_;
    /// Undo Button
    QPushButton * undobutton_;

    /// The tree.
    QTreeWidget * treeview_;
  };
}

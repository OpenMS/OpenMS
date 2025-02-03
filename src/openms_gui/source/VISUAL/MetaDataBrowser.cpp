// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


//OpenMS
#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/VISUAL/VISUALIZER/SampleVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>
#include <OpenMS/VISUAL/VISUALIZER/DigestionVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ModificationVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/TaggingVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/HPLCVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/GradientVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/SoftwareVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/SourceFileVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ContactPersonVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/InstrumentVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/IonSourceVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/IonDetectorVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/MassAnalyzerVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/DataProcessingVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ProteinIdentificationVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ProteinHitVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/PeptideHitVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ExperimentalSettingsVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/AcquisitionVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/AcquisitionInfoVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoDescriptionVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/PrecursorVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ProductVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/InstrumentSettingsVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/PeptideIdentificationVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/SpectrumSettingsVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/DocumentIdentifierVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ScanWindowVisualizer.h>

#include <QtWidgets/QLabel>
#include <QtWidgets/QStackedWidget>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QSplitter>

class QLayoutItem;

using namespace std;


//class QtWidgets/QLayoutItem;

namespace OpenMS
{

  MetaDataBrowser::MetaDataBrowser(bool editable, QWidget * parent, bool modal) :
    QDialog(parent),
    editable_(editable)
  {
    setWindowTitle("Meta data");

    setModal(modal);

    //splitter to hold treewidget and the right part (on a dummy widget)
    QGridLayout * grid = new QGridLayout(this);
    QSplitter * splitter = new QSplitter(Qt::Horizontal, this);
    grid->addWidget(splitter, 0, 0);

    //Create the tree for exploring data
    treeview_ = new QTreeWidget(this);
    treeview_->setColumnCount(3);
    treeview_->setHeaderLabel("Browse in Metadata tree");
    treeview_->setRootIsDecorated(true);
    treeview_->setColumnHidden(1, true);
    treeview_->setColumnHidden(2, true);
    splitter->addWidget(treeview_);

    //create a dummy widget to hold a grid layout and the item stack
    QWidget * dummy = new QWidget(splitter);
    splitter->addWidget(dummy);
    grid = new QGridLayout(dummy);
    grid->setColumnStretch(0, 1);

    //Create WidgetStack for managing visible metadata
    ws_ = new QStackedWidget(dummy);
    grid->addWidget(ws_, 0, 0, 1, 3);

    if (isEditable())
    {
      saveallbutton_ = new QPushButton("OK", dummy);
      cancelbutton_ = new QPushButton("Cancel", dummy);
      grid->addWidget(saveallbutton_, 1, 1);
      grid->addWidget(cancelbutton_, 1, 2);
      connect(saveallbutton_, SIGNAL(clicked()), this, SLOT(saveAll_()));
      connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject()));
    }
    else
    {
      closebutton_ = new QPushButton("Close", dummy);
      grid->addWidget(closebutton_, 1, 2);
      connect(closebutton_, SIGNAL(clicked()), this, SLOT(reject()));
    }

    connect(treeview_, SIGNAL(itemSelectionChanged()), this, SLOT(showDetails_()));

    status_list_ = "";

  }  //end of constructor

  void MetaDataBrowser::connectVisualizer_(BaseVisualizerGUI * ptr)
  {
    connect(ptr, SIGNAL(sendStatus(std::string)), this, SLOT(setStatus(std::string)));
  }

  bool MetaDataBrowser::isEditable() const
  {
    return editable_;
  }

  void MetaDataBrowser::setStatus(const std::string& status)
  {
    status_list_ = status_list_ + "\n" + status;
  }

  void MetaDataBrowser::showDetails_()
  {
    QList<QTreeWidgetItem *> list = treeview_->selectedItems();
    if (list.empty())
    {
      return;
    }
    ws_->setCurrentIndex(list[0]->text(1).toInt());
  }

  void MetaDataBrowser::saveAll_()
  {
    try
    {
      //call internal store function of all active visualizer objects
      for (int i = 0; i < ws_->count(); ++i)
      {
        dynamic_cast<BaseVisualizerGUI *>(ws_->widget(i))->store();
      }
      if (status_list_.length() != 0)
      {
        status_list_ = status_list_ + "\n" + "\n" + "Invalid modifications will not be saved.";
        QMessageBox::warning(this, tr("Save warning"), status_list_.c_str());
      }
    }
    catch (exception & e)
    {
      cout << "Exception while trying to save modifications." << endl << e.what() << endl;
    }

    //close dialog
    accept();
  }

  //-------------------------------------------------------------------------------
  // overloaded visualize functions to call the corresponding data visualizer
  // (in alphabetical oeder)
  //-------------------------------------------------------------------------------

  //Visualizing AcquisitionInfo object
  void MetaDataBrowser::visualize_(AcquisitionInfo & meta, QTreeWidgetItem * parent)
  {
    AcquisitionInfoVisualizer * visualizer = new AcquisitionInfoVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Acquisition Info" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    //Get Aquisition objects
    visualizeAll_(meta, item);

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing Acquisition object
  void MetaDataBrowser::visualize_(Acquisition & meta, QTreeWidgetItem * parent)
  {
    AcquisitionVisualizer * visualizer = new AcquisitionVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Acquisition" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing ContactPerson object
  void MetaDataBrowser::visualize_(ContactPerson & meta, QTreeWidgetItem * parent)
  {
    ContactPersonVisualizer * visualizer = new ContactPersonVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "ContactPerson" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing Digestion object
  void MetaDataBrowser::visualize_(Digestion & meta, QTreeWidgetItem * parent)
  {
    DigestionVisualizer * visualizer = new DigestionVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Digestion" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing ExperimentalSettings object
  void MetaDataBrowser::visualize_(ExperimentalSettings & meta, QTreeWidgetItem * parent)
  {
    ExperimentalSettingsVisualizer * visualizer = new ExperimentalSettingsVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "ExperimentalSettings" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<DocumentIdentifier &>(meta), item);

    //check for Sample
    visualize_(meta.getSample(), item);

    //check for ProteinIdentification
    visualizeAll_(meta.getProteinIdentifications(), item);

    //check for Instrument
    visualize_(meta.getInstrument(), item);

    //check for SourceFiles
    visualizeAll_(meta.getSourceFiles(), item);

    //check for ContactPersons
    visualizeAll_(meta.getContacts(), item);

    //check for HPLC
    visualize_(meta.getHPLC(), item);

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing Gradient object
  void MetaDataBrowser::visualize_(Gradient & meta, QTreeWidgetItem * parent)
  {
    GradientVisualizer * visualizer = new GradientVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Gradient" << QString::number(ws_->addWidget(visualizer));

    if (parent == nullptr)
    {
      new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      new QTreeWidgetItem(parent, labels);
    }
    connectVisualizer_(visualizer);
  }

  //Visualizing HPLC object
  void MetaDataBrowser::visualize_(HPLC & meta, QTreeWidgetItem * parent)
  {
    HPLCVisualizer * visualizer = new HPLCVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "HPLC" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }


    visualize_(meta.getGradient(), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing PeptideIdentification object
  void MetaDataBrowser::visualize_(PeptideIdentification & meta, QTreeWidgetItem * parent)
  {
    PeptideIdentificationVisualizer * visualizer = new PeptideIdentificationVisualizer(isEditable(), this, this);

    QStringList labels;
    int id = ws_->addWidget(visualizer);
    labels << QString("PeptideIdentification %1").arg(meta.getScoreType().c_str()) << QString::number(id);

    visualizer->load(meta, id);

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    //check for proteins and peptides hits
    meta.assignRanks();

    //list all peptides hits in the tree
    for (Size i = 0; i < meta.getHits().size(); ++i)
    {
      visualize_(const_cast<PeptideHit &>(meta.getHits()[i]), item);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing ProteinIdentification object
  void MetaDataBrowser::visualize_(ProteinIdentification & meta, QTreeWidgetItem * parent)
  {
    ProteinIdentificationVisualizer * visualizer = new ProteinIdentificationVisualizer(isEditable(), this, this);

    QStringList labels;
    int id = ws_->addWidget(visualizer);
    labels << QString("ProteinIdentification %1").arg(meta.getSearchEngine().c_str()) << QString::number(id);

    visualizer->load(meta, id);

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    //check for proteinhits objects
    meta.assignRanks();

    for (Size i = 0; i < meta.getHits().size(); ++i)
    {
      visualize_(const_cast<ProteinHit &>(meta.getHits()[i]), item);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing InstrumentSettings object
  void MetaDataBrowser::visualize_(InstrumentSettings & meta, QTreeWidgetItem * parent)
  {
    InstrumentSettingsVisualizer * visualizer = new InstrumentSettingsVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "InstrumentSettings" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    //ScanWindows
    visualizeAll_(meta.getScanWindows(), item);

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing Instrument object
  void MetaDataBrowser::visualize_(Instrument & meta, QTreeWidgetItem * parent)
  {
    InstrumentVisualizer * visualizer = new InstrumentVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Instrument" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    //IonSources
    visualizeAll_(meta.getIonSources(), item);

    //MassAnalyzers
    visualizeAll_(meta.getMassAnalyzers(), item);

    //IonDectector
    visualizeAll_(meta.getIonDetectors(), item);

    //Software
    visualize_(meta.getSoftware(), item);

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing IonDetector object
  void MetaDataBrowser::visualize_(IonDetector & meta, QTreeWidgetItem * parent)
  {
    IonDetectorVisualizer * visualizer = new IonDetectorVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "IonDetector" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing IonSource object
  void MetaDataBrowser::visualize_(IonSource & meta, QTreeWidgetItem * parent)
  {
    IonSourceVisualizer * visualizer = new IonSourceVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "IonSource" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing MassAnalyzer object
  void MetaDataBrowser::visualize_(MassAnalyzer & meta, QTreeWidgetItem * parent)
  {
    MassAnalyzerVisualizer * visualizer = new MassAnalyzerVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "MassAnalyzer" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing MetaInfoDescription object
  void MetaDataBrowser::visualize_(MetaInfoDescription & meta, QTreeWidgetItem * parent)
  {
    MetaInfoDescriptionVisualizer * visualizer = new MetaInfoDescriptionVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    String name = String("MetaInfoDescription ") + meta.getName();
    labels << name.c_str() << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    //check for DataProcessing
    visualizeAll_(meta.getDataProcessing(), item);

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing MetaInfoInterface object
  void MetaDataBrowser::visualize_(MetaInfoInterface & meta, QTreeWidgetItem * parent)
  {
    MetaInfoVisualizer * visualizer = new MetaInfoVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "MetaInfo" << QString::number(ws_->addWidget(visualizer));

    if (parent == nullptr)
    {
      new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      new QTreeWidgetItem(parent, labels);
    }
    connectVisualizer_(visualizer);
  }

  //Visualizing modification object
  void MetaDataBrowser::visualize_(Modification & meta, QTreeWidgetItem * parent)
  {
    ModificationVisualizer * visualizer = new ModificationVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Modification" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing PeptideHit object
  void MetaDataBrowser::visualize_(PeptideHit & meta, QTreeWidgetItem * parent)
  {
    PeptideHitVisualizer * visualizer = new PeptideHitVisualizer(isEditable(), this);
    visualizer->load(meta);

    String name = String("Pep ") + meta.getSequence().toString() + " (" + meta.getScore() + ')';
    QString qs_name(name.c_str());

    QStringList labels;
    labels << qs_name << QString::number(ws_->addWidget(visualizer)) << QString::number(meta.getScore());

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing Precursor object
  void MetaDataBrowser::visualize_(Precursor & meta, QTreeWidgetItem * parent)
  {
    PrecursorVisualizer * visualizer = new PrecursorVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Precursor" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing Product object
  void MetaDataBrowser::visualize_(Product & meta, QTreeWidgetItem * parent)
  {
    ProductVisualizer * visualizer = new ProductVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Product" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing DataProcessing object
  void MetaDataBrowser::visualize_(DataProcessingPtr & meta, QTreeWidgetItem * parent)
  {
    DataProcessingVisualizer * visualizer = new DataProcessingVisualizer(isEditable(), this);
    visualizer->load(*meta);

    QStringList labels;
    labels << "DataProcessing" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    //Software object
    visualize_(meta->getSoftware(), item);

    visualize_(dynamic_cast<MetaInfoInterface &>(*meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing ProteinHit object
  void MetaDataBrowser::visualize_(ProteinHit & meta, QTreeWidgetItem * parent)
  {
    ProteinHitVisualizer * visualizer = new ProteinHitVisualizer(isEditable(), this);
    visualizer->load(meta);

    String name = String("Prot ") + meta.getAccession() + " (" + meta.getScore() + ')';
    QString qs_name(name.c_str());

    QStringList labels;
    labels << qs_name << QString::number(ws_->addWidget(visualizer)) << QString::number(meta.getScore());

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing ProteinHit object
  void MetaDataBrowser::visualize_(ScanWindow & meta, QTreeWidgetItem * parent)
  {
    ScanWindowVisualizer * visualizer = new ScanWindowVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Scan window" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing sample object
  void MetaDataBrowser::visualize_(Sample & meta, QTreeWidgetItem * parent)
  {
    SampleVisualizer * visualizer = new SampleVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << (String("Sample ") + meta.getName()).c_str() << QString::number(ws_->addWidget(visualizer));
    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    //check for treatments
    if (meta.countTreatments() != 0)
    {
      for (Int i = 0; i < meta.countTreatments(); ++i)
      {
        if (meta.getTreatment(i).getType() == "Digestion")
        {
          visualize_((const_cast<Digestion &>(dynamic_cast<const Digestion &>(meta.getTreatment(i)))), item);
        }
        else if (meta.getTreatment(i).getType() == "Modification")
        {
          //Cast SampleTreatment reference to a const modification reference
          visualize_((const_cast<Modification &>(dynamic_cast<const Modification &>(meta.getTreatment(i)))), item);

        }
        else if (meta.getTreatment(i).getType() == "Tagging")
        {
          visualize_((const_cast<Tagging &>(dynamic_cast<const Tagging &>(meta.getTreatment(i)))), item);
        }
      }
    }

    //subsamples
    visualizeAll_(meta.getSubsamples(), item);

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing Software object
  void MetaDataBrowser::visualize_(Software & meta, QTreeWidgetItem * parent)
  {
    SoftwareVisualizer * visualizer = new SoftwareVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Software" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing SourceFile object
  void MetaDataBrowser::visualize_(SourceFile & meta, QTreeWidgetItem * parent)
  {
    SourceFileVisualizer * visualizer = new SourceFileVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "SourceFile" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }

    visualize_(dynamic_cast<MetaInfoInterface &>(meta), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing SpectrumSettings object
  void MetaDataBrowser::visualize_(SpectrumSettings & meta, QTreeWidgetItem * parent)
  {
    SpectrumSettingsVisualizer * visualizer = new SpectrumSettingsVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "SpectrumSettings" << QString::number(ws_->addWidget(visualizer));

    QTreeWidgetItem * item;
    if (parent == nullptr)
    {
      item = new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      item = new QTreeWidgetItem(parent, labels);
    }


    //check for InstrumentSettings
    visualize_(meta.getInstrumentSettings(), item);

    //check for DataProcessing
    visualizeAll_(meta.getDataProcessing(), item);

    //check for Precursors
    for (Size i = 0; i < meta.getPrecursors().size(); ++i)
    {
      visualize_(meta.getPrecursors()[i], item);
    }

    //check for Products
    for (Size i = 0; i < meta.getProducts().size(); ++i)
    {
      visualize_(meta.getProducts()[i], item);
    }

    //check for AcquisitionInfo
    visualize_(meta.getAcquisitionInfo(), item);

    //check for PeptideIdentification
    visualizeAll_(meta.getPeptideIdentifications(), item);

    connectVisualizer_(visualizer);
  }

  //Visualizing tagging object
  void MetaDataBrowser::visualize_(Tagging & meta, QTreeWidgetItem * parent)
  {
    TaggingVisualizer * visualizer = new TaggingVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "Tagging" << QString::number(ws_->addWidget(visualizer));

    if (parent == nullptr)
    {
      new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      new QTreeWidgetItem(parent, labels);
    }
    connectVisualizer_(visualizer);
  }

  //Visualizing tagging object
  void MetaDataBrowser::visualize_(DocumentIdentifier & meta, QTreeWidgetItem * parent)
  {
    DocumentIdentifierVisualizer * visualizer = new DocumentIdentifierVisualizer(isEditable(), this);
    visualizer->load(meta);

    QStringList labels;
    labels << "DocumentIdentifier" << QString::number(ws_->addWidget(visualizer));

    if (parent == nullptr)
    {
      new QTreeWidgetItem(treeview_, labels);
    }
    else
    {
      new QTreeWidgetItem(parent, labels);
    }
    connectVisualizer_(visualizer);
  }

  void MetaDataBrowser::filterHits_(double threshold, bool higher_better, int tree_item_id)
  {
    // find item in tree belonging to PeptideIdentification object
    QTreeWidgetItem * item = treeview_->findItems(QString::number(tree_item_id), Qt::MatchExactly | Qt::MatchRecursive, 1).first();

    //set the items visible or not visible depending to their score and the current threshold
    for (int i = 0; i < item->childCount(); ++i)
    {
      QTreeWidgetItem * child = item->child(i);

      if ((higher_better && child->text(2).toFloat() <= threshold) || (!higher_better && child->text(2).toFloat() >= threshold))
      {
        child->setHidden(true);
      }
      else
      {
        child->setHidden(false);
      }
    }

    //parent item must be collapsed and re-expanded so the items will be shown...
    treeview_->collapseItem(item);
    treeview_->expandItem(item);
  }

  /// Filters hits according to a score @a threshold. Takes the score orientation into account
  void MetaDataBrowser::showAllHits_(int tree_item_id)
  {
    // find item in tree belonging to PeptideIdentification object
    QTreeWidgetItem * item = treeview_->findItems(QString::number(tree_item_id), Qt::MatchExactly | Qt::MatchRecursive, 1).first();

    //set the items visible or not visible depending to their score and the current threshold
    for (int i = 0; i < item->childCount(); ++i)
    {
      item->child(i)->setHidden(false);
    }

    //parent item must be collapsed and re-expanded so the items will be shown...
    treeview_->collapseItem(item);
    treeview_->expandItem(item);
  }

  void MetaDataBrowser::add(Feature & feature)
  {
    //peptide ids
    for (PeptideIdentification& pep : feature.getPeptideIdentifications())
    {
      add(pep);
    }

    add(static_cast<MetaInfoInterface &>(feature));

    treeview_->expandItem(treeview_->findItems(QString::number(0), Qt::MatchExactly, 1).first());
  }

  void MetaDataBrowser::add(ConsensusFeature & feature)
  {
    //peptide ids
    for (PeptideIdentification& pep : feature.getPeptideIdentifications())
    {
      add(pep);
    }

    add(static_cast<MetaInfoInterface &>(feature));

    treeview_->expandItem(treeview_->findItems(QString::number(0), Qt::MatchExactly, 1).first());
  }

  void MetaDataBrowser::add(ConsensusMap & map)
  {
    //identifier
    add(static_cast<DocumentIdentifier &>(map));

    // protein identifications
    for (Size i = 0; i < map.getProteinIdentifications().size(); ++i)
    {
      add(map.getProteinIdentifications()[i]);
    }

    //unassigned peptide ids
    for (Size i = 0; i < map.getUnassignedPeptideIdentifications().size(); ++i)
    {
      add(map.getUnassignedPeptideIdentifications()[i]);
    }

    add(static_cast<MetaInfoInterface &>(map));

    treeview_->expandItem(treeview_->findItems(QString::number(0), Qt::MatchExactly, 1).first());
  }

} //end of MetaDataBrowser

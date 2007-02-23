// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPViewBase.h>

#ifdef DB_DEF
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/FORMAT/DB/DBAdapter.h>
#endif

#ifdef CGAL_DEF
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/VISUAL/DIALOGS/FeaFiDialog.h>
#endif

#include <OpenMS/VISUAL/DIALOGS/SaveImageDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewBasePDP.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWindow.h>
#include <OpenMS/VISUAL/Spectrum2DWindow.h>
#include <OpenMS/VISUAL/Spectrum3DWindow.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/VISUAL/DIALOGS/PeakPickingDialog.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolaySVDFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/VISUAL/DIALOGS/SmoothingDialog.h>
#include <OpenMS/FILTERING/BASELINE/TopHatFilter.h>
#include <OpenMS/VISUAL/DIALOGS/BaselineFilteringDialog.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>
#include <OpenMS/FORMAT/DFeaturePairsFile.h>

//Qt
#include <QtGui/QToolBar>
#include <QtGui/QDockWidget>
#include <QtGui/QListWidget>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPrinter>
#include <QtGui/QStatusBar>
#include <QtGui/QToolButton>
#include <QtGui/QMessageBox>
#include <QtGui/QListWidgetItem>
#include <QtGui/QToolTip>
#include <QtGui/QFileDialog>
#include <QtGui/QPainter>
#include <QtGui/QPrintDialog>
#include <QtCore/QDir>

//action modes
#include "../VISUAL/ICONS/zoom.xpm"
#include "../VISUAL/ICONS/translate.xpm"
#include "../VISUAL/ICONS/select.xpm"
#include "../VISUAL/ICONS/measure.xpm"

//intensity modes
#include "../VISUAL/ICONS/lin.xpm"
#include "../VISUAL/ICONS/log.xpm"
#include "../VISUAL/ICONS/percentage.xpm"
#include "../VISUAL/ICONS/snap.xpm"

//common
#include "../VISUAL/ICONS/grid.xpm"
#include "../VISUAL/ICONS/print.xpm"
#include "../VISUAL/ICONS/reset_zoom.xpm"
#include "../VISUAL/ICONS/tile_horizontal.xpm"
#include "../VISUAL/ICONS/tile_vertical.xpm"

//1d
#include "../VISUAL/ICONS/lines.xpm"
#include "../VISUAL/ICONS/peaks.xpm"

//2d
#include "../VISUAL/ICONS/points.xpm"
#include "../VISUAL/ICONS/colors.xpm"
#include "../VISUAL/ICONS/contours.xpm"

#include "../VISUAL/ICONS/TOPPView.xpm"

#include <algorithm>

using namespace std;

namespace OpenMS
{
  using namespace Internal;
	using namespace Math;

  TOPPViewBase::TOPPViewBase(QWidget* parent):
      QMainWindow(parent),
      PreferencesManager()
  {
  	setWindowTitle("TOPPView");
    setWindowIcon(QIcon(XPM_toppview));
    //prevents errors caused by too small width,height values
    setMinimumSize(400,400);

    // create dummy widget (to be able to have a layout), Tab bar and workspace
    QWidget* dummy = new QWidget(this);
    setCentralWidget(dummy);
    QVBoxLayout* box_layout = new QVBoxLayout(dummy);
    tab_bar_ = new EnhancedTabBar(dummy);
    tab_bar_->addTab("dummy");
    tab_bar_->setMinimumSize(tab_bar_->sizeHint());
    tab_bar_->removeTab(0);

    //connect slots and sigals for selecting spectra
    connect(tab_bar_,SIGNAL(selected(int)),this,SLOT(focusByTab(int)));
    connect(tab_bar_,SIGNAL(doubleClicked(int)),this,SLOT(closeByTab(int)));

    box_layout->addWidget(tab_bar_);
    ws_=new QWorkspace(dummy);
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateToolbar()));
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateTabBar(QWidget*)));
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateLayerbar()));
    box_layout->addWidget(ws_);

    // File menu
    QMenu* file = new QMenu("&File",this);
    menuBar()->addMenu(file);
    file->addAction("&Open",this,SLOT(openSpectrumDialog()));
    file->addAction("&Close",this,SLOT(closeFile()));
    file->addSeparator();
    
    QMenu* recent_menu = new QMenu("Recent files", this);
    //create the max mumber of recent files actions
  	recent_actions_.resize(20);
		for (Size i = 0; i<20; ++i)
		{
			recent_actions_[i] = recent_menu->addAction("",this,SLOT(openRecentFile()));
			recent_actions_[i]->setVisible(false);
		}
  	file->addMenu(recent_menu);

    file->addSeparator();
    file->addAction("&Preferences",this, SLOT(preferencesDialog()));
    file->addAction("&Quit",qApp,SLOT(quit()));
    
    //Layer menu
    QMenu* layer = new QMenu("&Layer",this);
    menuBar()->addMenu(layer);
    layer->addAction("&Save visible data",this,SLOT(saveLayer()));
    layer->addAction("&Edit metadata",this,SLOT(editMetadata()));
    layer->addAction("&Intensity distribution",this,SLOT(layerIntensityDistribution()));
		layer->addSeparator();
    layer->addAction("&Preferences",this, SLOT(layerPreferencesDialog()));
    
    //View menu
    QMenu * view = new QMenu("&View",this);
    menuBar()->addMenu(view);
    view->addAction("&Go to",this,SLOT(gotoDialog()), Qt::CTRL+Qt::Key_G);
   	view->addAction("Show/Hide axis legends",this,SLOT(changeAxisVisibility()));
   	
    //Image menu
    QMenu * image = new QMenu("&Image",this);
    menuBar()->addMenu(image);
    image->addAction("&Save to file",this,SLOT(saveImage()));
    image->addAction("&Print",this,SLOT(print()));

    //Windows menu
    QMenu * windows = new QMenu("&Windows", this);
    menuBar()->addMenu(windows);
    windows->addAction("&Cascade",this->ws_,SLOT(cascade()));
    windows->addAction("&Tile automatic",this->ws_,SLOT(tile()));
    windows->addAction(QIcon(QPixmap(XPM_tile_h)),"Tile &vertical",this,SLOT(tileHorizontal()));
    windows->addAction(QIcon(QPixmap(XPM_tile_v)),"Tile &horizontal",this,SLOT(tileVertical()));
    
    //Tools menu
    QMenu* tools_menu = new QMenu("&Tools",this);
    menuBar()->addMenu(tools_menu);
    tools_menu->addAction("&Show Selected Peaks (1D)", this, SLOT(showPeaklistActiveSpectrum()));
    tools_menu->addAction("&Pick Peaks", this, SLOT(pickActiveSpectrum()));
    tools_menu->addAction("&Smooth data", this, SLOT(smoothActiveSpectrum()));
    tools_menu->addAction("&Filter baseline in data", this, SLOT(baselineFilteringActiveSpectrum()));
    tools_menu->addAction("&Find Features (2D)", this, SLOT(findFeaturesActiveSpectrum()));

    //create status bar
    message_label_ = new QLabel(statusBar());
    statusBar()->addWidget(message_label_,1);

    rt_label_ = new QLabel("RT: 12345678", statusBar());
    rt_label_->setMinimumSize(rt_label_->sizeHint());
    rt_label_->setText("");
    statusBar()->addPermanentWidget(rt_label_,0);
    mz_label_ = new QLabel("m/z: 12345678", statusBar());
    mz_label_->setMinimumSize(mz_label_->sizeHint());
    mz_label_->setText("");
    statusBar()->addPermanentWidget(mz_label_,0);
    int_label_ = new QLabel("Int: 123456789012", statusBar());
    int_label_->setMinimumSize(int_label_->sizeHint());
    int_label_->setText("");
    statusBar()->addPermanentWidget(int_label_,0);

    //create toolbars and connect signals
    createToolBars_();

    //register meta value names
    MetaInfo::registry().registerName("FeatureDrawMode", "Specify what to draw of the Feature: BoundingBox, ConvexHulls");

    //set preferences file name + load preferencs
    loadPreferences();
  }

  TOPPViewBase::~TOPPViewBase()
  {
  	savePreferences();
  }

  void TOPPViewBase::closeEvent(QCloseEvent* /*event*/)
  {
  	QList<QWidget*> windows = ws_->windowList();
		for(int i=0; i< windows.size(); ++i)
		{
			windows.at(i)->close();
		}
  }

  void TOPPViewBase::addDBSpectrum(UnsignedInt db_id, bool as_new_window, bool maps_as_2d, bool maximize, OpenDialog::Mower use_mower)
  {
#ifdef DB_DEF
    //DBConnection for all DB queries
    DBConnection con;
    con.connect(getPref("Preferences:DB:Name"), getPref("Preferences:DB:Login"),getPref("DBPassword"),getPref("Preferences:DB:Host"),getPrefAsInt("Preferences:DB:Port"));

    //DB adapter
    DBAdapter dba(con);

    String db_id_string(db_id);
    QSqlQuery result;
    con.executeQuery("SELECT count(id) from DATA_Spectrum where fid_MSExperiment='"+db_id_string+"' and MSLevel='1'",result);

    //tab caption
    String caption = "DB ("+db_id_string+")";

    SpectrumWindow* w;

    SpectrumCanvas::ExperimentType* exp = 0;

    if (activeWindow_()==0)
    {
      as_new_window = true;
    }

    //open in new window
    if (as_new_window)
    {
      //create 1D View
      if (result.value(0).toInt()==1)
      {
      	//cout << "NEW 1D" << endl;
        // create 1D window
        w = new Spectrum1DWindow(ws_);

        //determine Spectrum id
        con.executeQuery("SELECT id from DATA_Spectrum where fid_MSExperiment='"+db_id_string+"' and MSLevel='1'",result);
        UID spectrum_id = result.value(0).toInt();

        //load data
        exp = &(w->widget()->canvas()->addEmptyPeakLayer());
        exp->resize(1);
        dba.loadSpectrum(spectrum_id,(*exp)[0]);
      }
      //create 2D/3D view
      else
      {
        //create 2D view
        if (maps_as_2d)
        {
          //cout << "NEW 2D" << endl;
          //create 2D window
          w = new Spectrum2DWindow(ws_);

          //load spectrum
          exp = &(w->widget()->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, (*exp));
        }
        //create 3D view
        else
        {
          //cout << "NEW 3D" << endl;
        	// create 3D window
          w = new Spectrum3DWindow(ws_);

          //load data
          exp = &(w->widget()->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, *exp);
        }
      }
      w->setMainPreferences(prefs_);
    }
    //open in active window
    else
    {
      //create 1D View
      if (result.value(0).toInt()==1)
      {
      	//cout << "ACTIVE 1D" << endl;
        w = active1DWindow_();
        //wrong active window type
        if (w==0)
        {
          QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+db_id_string+") in a 2D window!<BR>Please open the file in new tab.").c_str());
          return;
        }
        else //open it
        {
          //determine Spectrum id
          con.executeQuery("SELECT id from DATA_Spectrum where fid_MSExperiment='"+db_id_string+"' and MSLevel='1'", result);
          UID spectrum_id = result.value(0).toInt();

          //load data
	        exp = &(w->widget()->canvas()->addEmptyPeakLayer());
	        exp->resize(1);
	        dba.loadSpectrum(spectrum_id,(*exp)[0]);
        }
      }
      //create 2D/3D view
      else
      {
        Spectrum2DWindow* w2 = active2DWindow_();
        Spectrum3DWindow* w3 = active3DWindow_();
        //wrong active window type
        if (w2==0 && w3==0)
        {
          QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+db_id_string+") in a 2D/3D window!<BR>Please open the file in new tab.").c_str());
          return;
        }
        //create 2D view
        if (w2!=0)
        {
        	//cout << "ACTIVE 2D" << endl;
          w = w2;

          //load data
          exp = &(w->widget()->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, (*exp));
        }
        //create 3D view
        else
        {
        	//cout << "ACTIVE 3D" << endl;
          w = w3;

          //load data
          exp = &(w->widget()->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, (*exp));
        }
      }
		}
		
    //noise estimator
    float cutoff = 0;
    if(use_mower!=OpenDialog::NO_MOWER && exp->size()>1)
    {
      cutoff = estimateNoise_(*exp);
      cout << "Estimated noise level: " << cutoff << endl;
    }
    w->widget()->canvas()->finishAdding(cutoff);
		w->widget()->canvas()->setCurrentLayerName(caption);
    //use_mower

    //do for all windows
    if (as_new_window)
    {
      connectWindowSignals_(w);
      w->setWindowTitle(caption.c_str());
      addClient(w,caption);
      addTab_(w,caption);
    }

    //do for all (in active and in new window, 1D/2D/3D)
    if(maximize)
    {
      w->showMaximized();
    }

    //do for all windows
    updateLayerbar();
#endif
  }

  float TOPPViewBase::estimateNoise_(const SpectrumCanvas::ExperimentType& exp)
  {
    float noise = 0.0;
    UnsignedInt count = 0;
    srand(time(0));
    //cout << "size: " << exp.size() << endl;
    while (count<10)
    {
      UnsignedInt scan = (UnsignedInt)( (double)rand() / ((double)(RAND_MAX)+1.0f) * (double)(exp.size()-1) );

      if (scan < exp.size() && exp[scan].getMSLevel()==1)
      {
        vector<float> tmp;
        tmp.reserve(exp[scan].size());
        for(SpectrumCanvas::ExperimentType::SpectrumType::ConstIterator it = exp[scan].begin()
            ; it != exp[scan].end()
            ; ++it)
        {
          tmp.push_back(it->getIntensity());
        }
        std::sort(tmp.begin(),tmp.end());
        //cout << "scan: "<< scan <<" Groesse: " << tmp.size() << " Index: " << (UnsignedInt)ceil((float)(tmp.size()-1)/1.25f) << " Wert: "<< tmp[(UnsignedInt)ceil((float)(tmp.size()-1)/1.25f)] << endl;
        noise += tmp[(UnsignedInt)ceil((float)(tmp.size()-1)/1.25f)];

        ++count;
      }
    }
    noise = noise / 10.0f;
    return noise;
  }

  void TOPPViewBase::preferencesDialog()
  {
    setActive(true);
    if (showPreferencesDialog())
    {
      savePreferences();
    }
  }

  void TOPPViewBase::addSpectrum(const String& filename,bool as_new_window, bool maps_as_2d, bool maximize, OpenDialog::Mower use_mower, FileHandler::Type force_type)
  {
    if (!File::exists(filename))
    {
      QMessageBox::warning(this,"Open file error",("The file '"+filename+"' does not exist!").c_str());
      return;
    }

    //extract the filename without path
    String caption = File::basename(filename);

    //update recent files list
    addRecentFile_(filename);

    //windowpointer
    SpectrumWindow* w=0;
    
    if (activeWindow_()==0)
    {
      as_new_window = true;
    }

		//determine file type if not forced
		FileHandler fh;
		if (force_type==FileHandler::UNKNOWN)
		{
			if (force_type==FileHandler::UNKNOWN)
			{
				force_type = fh.getTypeByFileName(filename);
			}

			if (force_type==FileHandler::UNKNOWN)
			{
				force_type = fh.getTypeByContent(filename);
			}
		}

		if (force_type==FileHandler::UNKNOWN)
		{
      QMessageBox::warning(this,"Open file error",("Could not determine file type of '"+filename+"'!").c_str());
      return;
		}

    if (as_new_window)
    {
      if (force_type==FileHandler::DTA)
      {
        w = new Spectrum1DWindow(ws_);
      }
      else if (maps_as_2d || force_type==FileHandler::FEATURE || force_type==FileHandler::FEATURE_PAIRS) //2d or features
      {
        w = new Spectrum2DWindow(ws_);
      }
      else //3d
      {
        w = new Spectrum3DWindow(ws_);
      }

      //set main preferences
      w->setMainPreferences(prefs_);
    }
    else //!as_new_window
    {
      if (active1DWindow_()!=0) //1D window
      {
        if (force_type!=FileHandler::DTA)
        {
          QMessageBox::warning(this,"Wrong file type",("You cannot open 2D data ("+filename+") in a 1D window!").c_str());
          return;
        }

        w = active1DWindow_();
      }
      else if (active2DWindow_()!=0) //2d window
      {
        if (force_type==FileHandler::DTA)
        {
          QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+filename+") in a 2D window!").c_str());
          return;
        }
        w = active2DWindow_();
      }
      else if (active3DWindow_()!=0)//3d window
      {
        if (force_type==FileHandler::DTA)
        {
          QMessageBox::warning(this,"Wrong file type",("You cannot open 1D data ("+filename+") in a 3D window!").c_str());
          return;
        }
        w = active3DWindow_();
      }
    }

    //try to read the data from file
    if (force_type==FileHandler::FEATURE) //features
    {
      DFeatureMap<2> map;
      try
      {
        DFeatureMapFile().load(filename,map);
      }
      catch(Exception::Base& e)
      {
        QMessageBox::warning(this,"Error",(String("Error while reading feature file: ")+e.what()).c_str());
        return;
      }
      w->widget()->canvas()->addLayer(map,false);
      w->widget()->canvas()->setCurrentLayerName(caption);
    }
    else if (force_type==FileHandler::FEATURE_PAIRS) //feature pairs
    {
    	//load pairs
      DFeaturePairVector<2> pairs;
      try
      {
        DFeaturePairsFile().load(filename,pairs);
      }
      catch(Exception::Base& e)
      {
        QMessageBox::warning(this,"Error",(String("Error while reading feature pairs file: ")+e.what()).c_str());
        return;
      }
      
      //convert to features
      DFeatureMap<2> map;
      DFeaturePairsFile::pairsToFeatures(pairs,map);
      w->widget()->canvas()->addLayer(map,true);
      w->widget()->canvas()->setCurrentLayerName(caption);
    }
    else
    {
      //try to read the data from file (raw/peak data)
      SpectrumCanvas::ExperimentType* exp = &(w->widget()->canvas()->addEmptyPeakLayer());
      try
      {
        fh.loadExperiment(filename,*exp, force_type);
      }
      catch(Exception::Base& e)
      {
        QMessageBox::warning(this,"Error",(String("Error while reading data file: ")+e.what()).c_str());
        return;
      }

      //check if only one scan is in a 2d file
      if (as_new_window && active1DWindow_()==0 && exp->size()==1)
      {
        delete(w);
        w = new Spectrum1DWindow(ws_);
        w->setMainPreferences(prefs_);
        exp = &(w->widget()->canvas()->addEmptyPeakLayer());
        FileHandler().loadExperiment(filename,*exp, force_type);
      }

      //do for all (in active and in new window, 1D/2D/3D)
      float cutoff=0;

      //calculate noise
      if(use_mower!=OpenDialog::NO_MOWER && exp->size()>1)
      {
        cutoff = estimateNoise_(*exp);
      }
      w->widget()->canvas()->finishAdding(cutoff);
      w->widget()->canvas()->setCurrentLayerName(caption);
    }
    
  	updateLayerbar();

    if (as_new_window)
    {
    	ws_->addWindow(w);
      connectWindowSignals_(w);
      w->setWindowTitle(filename.c_str());
      addClient(w,filename);
      addTab_(w,caption);
    }
    if(maximize)
    {
      w->showMaximized();
    }
  }

  void TOPPViewBase::addRecentFile_(const String& filename)
  {
    String tmp = filename;

    //add prefix to relative paths
    if (!tmp.hasPrefix("/"))
    {
      tmp = QDir::currentPath().toAscii().data()+string("/")+ tmp;
    }
	
		// remove the new file if already in the recent list and prepend it
		recent_files_.removeAll(tmp.c_str());
		recent_files_.prepend(tmp.c_str());
		
		//remove those files exceeding the defined number
		UnsignedInt number_of_recent_files = UnsignedInt(prefs_.getValue("Preferences:NumberOfRecentFiles"));
		while ((Size)recent_files_.size() > number_of_recent_files)
		{
			recent_files_.removeLast();
		}
		
    updateRecentMenu_();
  }

  void TOPPViewBase::updateRecentMenu_()
  {
    //get/correct number of recent files
		UnsignedInt number_of_recent_files = UnsignedInt(prefs_.getValue("Preferences:NumberOfRecentFiles"));
		if (number_of_recent_files>20)
		{
			number_of_recent_files = 20;
			prefs_.setValue("Preferences:NumberOfRecentFiles",20);
		}
		
		for (UnsignedInt i = 0; i < 20; ++i)
		{
			if (i < (UnsignedInt)(recent_files_.size()))
			{
				QString text = tr("&%1 %2").arg(i).arg(recent_files_[i]);
				recent_actions_[i]->setText(text);
				recent_actions_[i]->setData(recent_files_[i]);
				recent_actions_[i]->setVisible(true);
			}
			else
			{
				recent_actions_[i]->setVisible(false);
			}
		}
  }

  void TOPPViewBase::maximizeActiveSpectrum()
  {
    if (ws_->activeWindow())
    {
      ws_->activeWindow()->showMaximized();
    }
  }

  void TOPPViewBase::addTab_(SpectrumWindow* w, const String& tabCaption)
  {
  	//static window counter
  	static int window_counter = 4711;
  	w->window_id = ++window_counter;

		//add tab and assign window id as data  	
    int tab_index = tab_bar_->addTab(tabCaption.c_str());
    tab_bar_->setTabData(tab_index, window_counter);
    //cout << "Added tab: '" << tabCaption << "' => " << window_counter << endl;
        
    //connect slots and sigals for removing the spectrum from the bar, when it is closed
    connect(w,SIGNAL(aboutToBeDestroyed(int)),this,SLOT(removeTab(int)));
    tab_bar_->setCurrentIndex(tab_index);
  }

  SpectrumWindow* TOPPViewBase::window_(int id) const
  {
  	//cout << "Looking for tab with id: " << id << endl;
  	QList<QWidget*> windows = ws_->windowList();
		for(int i=0; i< windows.size(); ++i)
		{
			SpectrumWindow* window = dynamic_cast<SpectrumWindow*>(windows.at(i));
			//cout << "  Tab " << i << ": " << window->window_id << endl;
			if (window->window_id == id)
			{
				return window;
			}
		}
		return 0;
  }

  void TOPPViewBase::removeTab(int id)
  {
  	for (int i=0; i < tab_bar_->count(); ++i)
  	{ 
  		if( tab_bar_->tabData(i).toInt() == id)
  		{
  			tab_bar_->removeTab(i);
  			return;
  		}
  	}
  }

  void TOPPViewBase::closeByTab(int index)
  {
  	SpectrumWindow* window = window_(tab_bar_->tabData(index).toInt());
  	if (window)
  	{
  		window->close();
  	}
  }
 
  void TOPPViewBase::focusByTab(int index)
  {
  	SpectrumWindow* window = window_(tab_bar_->tabData(index).toInt());
  	if (window)
  	{
  		window->setFocus();
  	}
  }

  void TOPPViewBase::saveImage()
  {
    //check if there is a active window
    SpectrumWindow* window = activeWindow_();
    if (window!=0)
    {
      SaveImageDialog* dialog = new SaveImageDialog(this);
      dialog->setSize(1024,768);
      if (dialog->exec())
      {
        QString format=dialog->getFormat();
        QString file_name = QFileDialog::getSaveFileName(this, "Save file", window->windowTitle().section('.',0,0),
                            format+" Images(*."+format+" *."+format.toLower()+")");
        if (!file_name.isEmpty())
        {
          //append missing file extension
          if (!file_name.toUpper().endsWith(format))
          {
            file_name.append("."+format.toLower());
          }

          QImage image = window->widget()->getImage(dialog->getXSize(),dialog->getYSize());
          image.save(file_name,format.toAscii().data(),100);
        }
      }
    }
  }

  void TOPPViewBase::print()
  {
#ifndef QT_NO_PRINTER
    SpectrumWindow* window = activeWindow_();
    if (window!=0)
    {
      QPrinter* printer = new QPrinter(QPrinter::HighResolution);
      printer->setResolution(300);
      
      QPrintDialog dialog(printer, this);
 			if (dialog.exec())
      {
        QPainter p;
        if (!p.begin(printer)) return;
        unsigned int dpix = printer->logicalDpiX();
        unsigned int dpiy = printer->logicalDpiY();
        QRect body(dpix,dpiy,printer->width()-2*dpix,printer->height()-2*dpiy);			// one inch margin
        QImage image = window->widget()->getImage(body.width(),body.height());
        p.drawImage(body,image);
        QString titel = QString("%1\n%2").arg(window->windowTitle().section('/',-1)).arg(QDate::currentDate().toString());
        //	p.drawText(dpix,0,body.width(),dpiy, Qt::AlignCenter, titel);
        p.drawText(dpix,body.height()+dpiy, body.width(), dpiy, Qt::AlignCenter, titel);
      }
      delete(printer);
    }
#endif

  }

  void TOPPViewBase::closeFile()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      ws_->activeWindow()->close();
    }
  }

  void TOPPViewBase::saveLayer()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      const LayerData& layer = activeWindow_()->widget()->canvas()->getCurrentLayer();
      //warn if hidden layer => wrong layer selected...
    	if (!layer.visible)
    	{
    		QMessageBox::warning(this,"Warning","The current layer is not visible!");
    	}
    	//Visible area
    	const SpectrumCanvas::AreaType& area = activeWindow_()->widget()->canvas()->getVisibleArea();
    	
    	if (layer.type==LayerData::DT_PEAK)
    	{
    		//Extract selected visible data to out
    		LayerData::ExperimentType out;
    		out.ExperimentalSettings::operator=(layer.peaks);
    		LayerData::ExperimentType::ConstIterator begin;
    		LayerData::ExperimentType::ConstIterator end; 
    		if (layer.peaks.size()==1)
    		{
	    		begin = layer.peaks.begin();
	    		end = layer.peaks.end();
    		}
    		else
    		{
	    		begin = layer.peaks.RTBegin(area.min()[1]);
	    		end = layer.peaks.RTBegin(area.max()[1]);
    		}
    		out.resize(end-begin);
				
				UnsignedInt i = 0;
    		for (LayerData::ExperimentType::ConstIterator it=begin; it!=end; ++it)
    		{
  				out[i].SpectrumSettings::operator=(*it);
  				out[i].setRetentionTime(it->getRetentionTime());
  				out[i].setMSLevel(it->getMSLevel());
  				out[i].setPrecursorPeak(it->getPrecursorPeak());
  				for (LayerData::ExperimentType::SpectrumType::ConstIterator it2 = it->MZBegin(area.min()[0]); it2!= it->MZEnd(area.max()[0]); ++it2)
  				{
  					if ( it2->getIntensity() >= layer.min_int && it2->getIntensity() <= layer.max_int)
  					{
  						out[i].push_back(*it2);
  					}
  				}
  				++i;
    		}
    		//no extracted data
    		if (out.size()==0)
    		{
    		  QMessageBox::warning(this,"Warning","The displayed region of the current layer is empty!");
    		  return;
    		}
    		//one scan => DTA
    		else if (out.size()==1)
    		{
		      QString file_name = QFileDialog::getSaveFileName(this, "Save file", prefs_.getValue("Preferences:DefaultPath").toString().c_str(),
					                    "DTA files (*.dta)");
					if (!file_name.isEmpty())
					{
					  DTAFile().store(file_name.toAscii().data(),out[0]);
					}
    		}
    		//more than one scan => MzData
    		else
    		{
		      QString file_name = QFileDialog::getSaveFileName(this, "Save file",  prefs_.getValue("Preferences:DefaultPath").toString().c_str(),
					                    "MzData files (*.mzData)");
					if (!file_name.isEmpty())
					{
					  MzDataFile().store(file_name.toAscii().data(),out);
					}
    		}
			}
			else //feature data
			{
		    enum DimensionId
		    {
	        RT = DimensionDescription < LCMS_Tag >::RT,
	        MZ = DimensionDescription < LCMS_Tag >::MZ
		    };
	    
    		//Extract selected visible data to out
    		LayerData::FeatureMapType out;
    		out.ExperimentalSettings::operator=(layer.features);
    		for (LayerData::FeatureMapType::ConstIterator it=layer.features.begin(); it!=layer.features.end(); ++it)
    		{
					if ( it->getIntensity() >= layer.min_int && 
							 it->getIntensity() <= layer.max_int &&
							 it->getPosition()[RT] >= area.min()[1] &&
							 it->getPosition()[RT] <= area.max()[1] &&
							 it->getPosition()[MZ] >= area.min()[0] &&
							 it->getPosition()[MZ] <= area.max()[0]
						 )
					{
						out.push_back(*it);
					}
  			}
    		//no extracted data
    		if (out.size()==0)
    		{
    		  QMessageBox::warning(this,"Warning","The displayed region of the current layer is empty!");
    		  return;
    		}
    		else
    		{
		      QString file_name = QFileDialog::getSaveFileName(this, "Save file",  prefs_.getValue("Preferences:DefaultPath").toString().c_str(),
					                    "features files (*.feat)" );
					if (!file_name.isEmpty())
					{
					  DFeatureMapFile().store(file_name.toAscii().data(),out);
					}
    		}
			}
    }
  }

  void TOPPViewBase::editMetadata()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      const LayerData& layer = activeWindow_()->widget()->canvas()->getCurrentLayer();
      //warn if hidden layer => wrong layer selected...
    	if (!layer.visible)
    	{
    		QMessageBox::warning(this,"Warning","The current layer is not visible!");
    	}
			MSMetaDataExplorer dlg(true, this);
      dlg.setWindowTitle("Edit meta data");
			if (layer.type==LayerData::DT_PEAK) //peak data
    	{
    		dlg.add(&(const_cast<LayerData&>(layer).peaks));
    	}
    	else //feature data
    	{
    		dlg.add(&(const_cast<LayerData&>(layer).features));
    	}
      dlg.exec();
    }
  }

  void TOPPViewBase::layerPreferencesDialog()
  {
    if (ws_->activeWindow())
    {
      activeWindow_()->setActive(true);
	  	if (showPreferencesDialog())
	    {
	      savePreferences();
	    }
    }
  }

  void TOPPViewBase::layerIntensityDistribution()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      activeWindow_()->widget()->showIntensityDistribution();
    }  	
  }

  void TOPPViewBase::changeAxisVisibility()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      activeWindow_()->widget()->showLegend(!activeWindow_()->widget()->isLegendShown());
    }  	
  }

  void TOPPViewBase::createToolBars_()
  {
  	QToolButton* b;
  	
  	//**Basic tool bar for all views**
    tool_bar_ = addToolBar("Basic tool bar");
    
    //action modes
    action_group_ = new QButtonGroup(tool_bar_);
    action_group_->setExclusive(true);
    
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(XPM_zoom));
    b->setToolTip("Action: Zoom");
    b->setShortcut(Qt::Key_Z);
    b->setCheckable(true);
    action_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::AM_ZOOM);
		tool_bar_->addWidget(b);
		
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(XPM_translate));
    b->setToolTip("Action: Translate");
    b->setShortcut(Qt::Key_T);
    b->setCheckable(true);
    action_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::AM_TRANSLATE);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(XPM_select));
    b->setToolTip("Action: Select");
    b->setShortcut(Qt::Key_S);
    b->setCheckable(true);
    action_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::AM_SELECT);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(XPM_measure));
    b->setToolTip("Action: Measure");
    b->setShortcut(Qt::Key_M);
    b->setCheckable(true);
    action_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::AM_MEASURE);
		tool_bar_->addWidget(b);

    connect(action_group_,SIGNAL(buttonClicked(int)),this,SLOT(setActionMode(int)));
    tool_bar_->addSeparator();
	   
    //intensity modes
    intensity_group_ = new QButtonGroup(tool_bar_);
    intensity_group_->setExclusive(true);
    
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(XPM_lin));
    b->setToolTip("Intensity: Normal");
    b->setShortcut(Qt::Key_N);
    b->setCheckable(true);
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_NONE);
		tool_bar_->addWidget(b);
    
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(XPM_log));
    b->setToolTip("Intensity: Logarithmic");
    b->setShortcut(Qt::Key_L);
    b->setCheckable(true);
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_LOG);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(XPM_percentage));
    b->setToolTip("Intensity: Percentage");
    b->setShortcut(Qt::Key_P);
    b->setCheckable(true);
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_PERCENTAGE);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(XPM_snap));
    b->setToolTip("Intensity: Snap to maximum displayed intensity");
    b->setShortcut(Qt::Key_A);
    b->setCheckable(true);
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_SNAP);
		tool_bar_->addWidget(b);

    connect(intensity_group_,SIGNAL(buttonClicked(int)),this,SLOT(setIntensityMode(int)));
    tool_bar_->addSeparator();

    //common buttons
    QAction* reset_zoom_button = tool_bar_->addAction(QPixmap(XPM_reset_zoom), "Reset Zoom", this, SLOT(resetZoom()));
    reset_zoom_button->setShortcut(Qt::Key_Backspace);
    tool_bar_->addSeparator();

    grid_button_ = tool_bar_->addAction(QPixmap(XPM_grid), "Show grid");
    grid_button_->setIcon(QIcon(QPixmap(XPM_grid)));
    grid_button_->setCheckable(true);
    connect(grid_button_,SIGNAL(toggled(bool)),this,SLOT(showGridLines(bool)));
    tool_bar_->addSeparator();

    QAction* print_button = tool_bar_->addAction(QPixmap(XPM_print), "Print", this, SLOT(print()));
    print_button->setShortcut(Qt::CTRL + Qt::Key_P);

    tool_bar_->resize(tool_bar_->sizeHint());
    tool_bar_->show();

    //**1D toolbar**
    tool_bar_1d_ = addToolBar("1D tool bar");

    //draw modes 1D
    draw_group_1d_ = new QButtonGroup(tool_bar_1d_);
    draw_group_1d_->setExclusive(true);
    
    b = new QToolButton(tool_bar_1d_);
    b->setIcon(QPixmap(XPM_peaks));
    b->setToolTip("Peak mode");
    b->setShortcut(Qt::Key_I);
    b->setCheckable(true);
    draw_group_1d_->addButton(b,Spectrum1DCanvas::DM_PEAKS);
		tool_bar_1d_->addWidget(b);
    
    b = new QToolButton(tool_bar_1d_);
    b->setIcon(QPixmap(XPM_lines));
    b->setToolTip("Raw data mode");
    b->setShortcut(Qt::Key_R);
    b->setCheckable(true);
    draw_group_1d_->addButton(b,Spectrum1DCanvas::DM_CONNECTEDLINES);
		tool_bar_1d_->addWidget(b);

    connect(draw_group_1d_,SIGNAL(buttonClicked(int)),this,SLOT(setDrawMode1D(int)));
    tool_bar_->addSeparator();


    link_box_ = new QComboBox(tool_bar_1d_);
    link_box_->setToolTip("Use this combobox to link two spectra.<BR>Linked spectra zoom in/out together.");
    tool_bar_1d_->addWidget(link_box_);
    connect(link_box_,SIGNAL(activated(int)),this,SLOT(linkActiveTo(int)));

    //**2D toolbar**
    tool_bar_2d_ = addToolBar("2D tool bar");

    dm_points_2d_ = tool_bar_2d_->addAction(QPixmap(XPM_points),"Show dots");
    dm_points_2d_->setShortcut(Qt::Key_D);
    dm_points_2d_->setCheckable(true);
    connect(dm_points_2d_, SIGNAL(toggled(bool)), this, SLOT(showPoints(bool)));

    dm_surface_2d_ = tool_bar_2d_->addAction(QPixmap(XPM_colors),"Show colored surface");
    dm_surface_2d_->setShortcut(Qt::Key_U);
    dm_surface_2d_->setCheckable(true);
    connect(dm_surface_2d_, SIGNAL(toggled(bool)), this, SLOT(showSurface(bool)));

    dm_contours_2d_ = tool_bar_2d_->addAction(QPixmap(XPM_contours),"Show contour lines");
    dm_contours_2d_->setShortcut(Qt::Key_C);
    dm_contours_2d_->setCheckable(true);
    connect(dm_contours_2d_, SIGNAL(toggled(bool)), this, SLOT(showContours(bool)));

    //**layer bar**
    QDockWidget* layer_bar = new QDockWidget("Layers", this);
    addDockWidget(Qt::RightDockWidgetArea, layer_bar);
    layer_manager_ = new QListWidget(layer_bar);
    layer_bar->setWidget(layer_manager_);
    layer_manager_->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(layer_manager_,SIGNAL(currentRowChanged(int)),this,SLOT(layerSelectionChange(int)));
		connect(layer_manager_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(layerContextMenu(const QPoint&)));
		connect(layer_manager_,SIGNAL(itemChanged(QListWidgetItem*)),this,SLOT(layerVisibilityChange(QListWidgetItem*)));
  }

  void TOPPViewBase::linkActiveTo(int index)
  {
  	Spectrum1DWindow* active = active1DWindow_();
  	if (active==0) return;
  	
  	//cout << "linkActiveTo() active: " << active->window_id << endl;
  	
  	//remove link if present
  	if (link_map_.find(active->window_id)!=link_map_.end())
  	{
  		SpectrumWindow* linked_to = window_(link_map_[active->window_id]);
  		if (linked_to != 0)
  		{
  			//cout << "  disconnect signals to: " << linked_to->window_id << endl;
	  		//disconnect outgoing signals
	  		disconnect(active->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),linked_to->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
  			//disconnect incoming signals
  			disconnect(linked_to->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),active->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
  			//erase entries in id map
  			link_map_.erase(active->window_id);
  			link_map_.erase(linked_to->window_id);
  		}
  		else
  		{
  			cout << "linkActiveTo() - disconnect: Error, could not find window with id '" << active->window_id << "'!" << endl;
  		}
  	}
  	  
    if (link_box_->itemText(index)=="<unlinked>") return;
    
    //link
    SpectrumWindow* link_to = window_(link_box_->itemData(index).toInt());
		if (link_to!=0)
		{
			//cout << "  connecting signals to: " << link_to->window_id << endl;
		  connect(active->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),link_to->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
      connect(link_to->widget()->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),active->widget()->canvas(),SLOT(setVisibleArea(DRange<2>)));
      //add entries to id map
      link_map_[link_to->window_id]=active->window_id;
      link_map_[active->window_id]=link_to->window_id;
		}
		else
		{
			cout << "linkActiveTo() - connect: Error, could not find window with id '" << link_box_->itemData(index).toInt() << "'!" << endl;
		}    
  }

  void TOPPViewBase::showStatusMessage(string msg, OpenMS::UnsignedInt time)
  {
    if (time==0)
    {
      message_label_->setText(msg.c_str());
      statusBar()->repaint();
    }
    else
    {
      statusBar()->showMessage(msg.c_str(), time);
    }
  }

  void TOPPViewBase::showCursorStatus(double mz, double intensity, double rt)
  {
    message_label_->setText("");
    if (mz==-1)
    {
      mz_label_->setText("m/z: ");
    }
    else
    {
      mz_label_->setText(("m/z: "+String(mz,8).fillLeft(' ',8)).c_str());
    }
    if (rt==-1)
    {
      rt_label_->setText("RT: ");
    }
    else
    {
      rt_label_->setText(("RT: "+String(rt,8).fillLeft(' ',8)).c_str());
    }
    if (intensity==-1)
    {
      int_label_->setText("Int: ");
    }
    else
    {
      int_label_->setText(("Int: "+String(intensity,12).fillLeft(' ',12)).c_str());
    }
    statusBar()->repaint();
  }

  void TOPPViewBase::showGridLines(bool b)
  {
    SpectrumWindow* window = activeWindow_();
    if (window!=0)
    {
      window->widget()->canvas()->showGridLines(b);
    }
  }

  void TOPPViewBase::resetZoom()
  {
    SpectrumWindow* window = activeWindow_();
    if (window!=0)
    {
      window->widget()->canvas()->resetZoom();
    }
  }

  void TOPPViewBase::setActionMode(int index)
  {
    SpectrumWindow* w = activeWindow_();
    if (w)
    {
    	w->widget()->setActionMode((OpenMS::SpectrumCanvas::ActionModes)index);
  	}
  }

  void TOPPViewBase::setIntensityMode(int index)
  {
    SpectrumWindow* w = activeWindow_();
    if (w)
    {
    	w->widget()->setIntensityMode((OpenMS::SpectrumCanvas::IntensityModes)index);
  	}
  }

  void TOPPViewBase::setDrawMode1D(int index)
  {
    Spectrum1DWindow* w = active1DWindow_();
    if (w)
    {
    	w->widget()->canvas()->setDrawMode((OpenMS::Spectrum1DCanvas::DrawModes)index);
  	}
  }

  void TOPPViewBase::showPoints(bool on)
  {
    if (Spectrum2DWindow* win = active2DWindow_())
    {
      win->widget()->canvas()->showPoints(on);
    }
  }

  void TOPPViewBase::showSurface(bool on)
  {
    if (Spectrum2DWindow* win = active2DWindow_())
    {
      win->widget()->canvas()->showSurface(on);
    }
  }

  void TOPPViewBase::showContours(bool on)
  {
    if (Spectrum2DWindow* win = active2DWindow_())
    {
      win->widget()->canvas()->showContours(on);
    }
  }

  void TOPPViewBase::updateToolbar()
  {
    SpectrumWindow* w = activeWindow_();

    if (w)
    {
      //set action mode
     	action_group_->button(w->widget()->getActionMode())->setChecked(true);

      //set intensity mode
     	intensity_group_->button(w->widget()->canvas()->getIntensityMode())->setChecked(true);

      //grid lines
      grid_button_->setChecked(w->widget()->canvas()->gridLinesShown());
    }

    //1D
    Spectrum1DWindow* w1 = active1DWindow_();
    if (w1)
    {
      //draw mode
      draw_group_1d_->button(w1->widget()->canvas()->getDrawMode())->setChecked(true);

      //update link selector
      int item_index = -1;
      link_box_->clear();
      link_box_->insertItem(++item_index,"<unlinked>",0);
      QWidgetList windows = ws_->windowList();
      for ( int i = 0; i < windows.count(); ++i )
      {
        Spectrum1DWindow* window = dynamic_cast<Spectrum1DWindow*>(windows.at(i));
        if (window !=0 && window!=w)
        {
          link_box_->insertItem(++item_index,File::basename(window->windowTitle().toAscii().data()).c_str(),window->window_id);
        	if (link_map_.find(w1->window_id)!=link_map_.end() && link_map_[w1->window_id] == window->window_id)
          {
            link_box_->setCurrentIndex(item_index);
          }
        }
      }

      //show/hide toolbars and buttons
      tool_bar_1d_->show();
      tool_bar_2d_->hide();
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_MEASURE)->setEnabled(false);
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_SELECT)->setEnabled(true);
    }

    //2d
    Spectrum2DWindow* w2 = active2DWindow_();
    if (w2)
    {
      //draw modes
      dm_points_2d_->setChecked(w2->widget()->canvas()->dotsAreShown());
      dm_surface_2d_->setChecked(w2->widget()->canvas()->surfaceIsShown());
      dm_contours_2d_->setChecked(w2->widget()->canvas()->contoursAreShown());

      //show/hide toolbars and buttons
      tool_bar_1d_->hide();
      tool_bar_2d_->show();
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_MEASURE)->setEnabled(true);
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_SELECT)->setEnabled(true);
    }

    //1D
    Spectrum3DWindow* w3 = active3DWindow_();
    if (w3)
    {
      //show/hide toolbars and buttons
      tool_bar_1d_->hide();
      tool_bar_2d_->hide();
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_MEASURE)->setEnabled(false);
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_SELECT)->setEnabled(false);
    }
  }

  void TOPPViewBase::updateLayerbar()
  {
		layer_manager_->clear();

    SpectrumCanvas* cc = activeCanvas_();
    if (cc == 0)
    {
      return;
    }
		QListWidgetItem* item = 0;
    for (UnsignedInt i = 0; i<cc->getLayerCount(); ++i)
    {
    	//add item
    	item = new QListWidgetItem( layer_manager_ );
			item->setText(cc->getLayer(i).name.c_str());
    	if (cc->getLayer(i).visible)
    	{
    		item->setCheckState(Qt::Unchecked);
    	}
    	else
    	{
    		item->setCheckState(Qt::Checked);
    	}
    	//highlight active item
    	if (i == cc->activeLayerIndex())
    	{
		    layer_manager_->blockSignals(true);
				layer_manager_->setCurrentItem(item);
				layer_manager_->blockSignals(false);
    	}
    }
  }

	void TOPPViewBase::layerSelectionChange(int i)
	{
		if (i!=-1) activeCanvas_()->activateLayer(i);
	}

	void TOPPViewBase::layerContextMenu(const QPoint & pos)
	{
		QListWidgetItem* item = layer_manager_->itemAt(pos);
		if (item && item!=layer_manager_->item(0))
		{
			QMenu* context_menu = new QMenu(layer_manager_);
			QAction* delete_action = context_menu->addAction("Delete");
			if (context_menu->exec(layer_manager_->mapToGlobal(pos)) == delete_action)
			{
				activeCanvas_()->removeLayer(layer_manager_->row(item));
				updateLayerbar();
			}
			delete (context_menu);
		}
	}

	void TOPPViewBase::layerVisibilityChange(QListWidgetItem* item)
	{
		int layer = layer_manager_->row(item);
		bool visible = activeCanvas_()->getLayer(layer).visible;
		
		if (item->checkState()==Qt::Unchecked && visible)
		{
			activeCanvas_()->changeVisibility(layer,false);
		}
		else if (item->checkState()==Qt::Checked && !visible)
		{
			activeCanvas_()->changeVisibility(layer,true);
		}
			
			
	}

  void TOPPViewBase::updateTabBar(QWidget* w)
  {
  	if (w)
  	{
  		tab_bar_->setCurrentIndex(dynamic_cast<SpectrumWindow*>(w)->window_id);
  	}
  }

  void TOPPViewBase::tileVertical()
  {
    // primitive horizontal tiling
    QWidgetList windows = ws_->windowList();
    if ( !windows.count() )
      return;

    if (active1DWindow_())
      active1DWindow_()->showNormal();
    if (active2DWindow_())
      active2DWindow_()->showNormal();

    int heightForEach = ws_->height() / windows.count();
    int y = 0;
    for ( int i = 0; i < int(windows.count()); ++i )
    {
      QWidget *window = windows.at(i);
      if ( window->isMaximized() || window->isFullScreen() )
      {

        // prevent flicker
        window->hide();
        window->setWindowState(Qt::WindowNoState);
        window->show();
      }
      int preferredHeight = window->minimumHeight()+window->parentWidget()->baseSize().height();
      int actHeight = std::max(heightForEach, preferredHeight);

      window->parentWidget()->setGeometry( 0, y, ws_->width(), actHeight );
      y += actHeight;
    }
  }

  void TOPPViewBase::tileHorizontal()
  {
    // primitive horizontal tiling
    QWidgetList windows = ws_->windowList();
    if ( !windows.count() )
      return;

    if (active1DWindow_())
      active1DWindow_()->showNormal();
    if (active2DWindow_())
      active2DWindow_()->showNormal();

    int widthForEach = ws_->width() / windows.count();
    int y = 0;
    for ( int i = 0; i < int(windows.count()); ++i )
    {
      QWidget *window = windows.at(i);
      if ( window->windowState() & Qt::WindowMaximized )
      {
        // prevent flicker
        window->hide();
        window->showNormal();
      }
      int preferredWidth = window->minimumWidth()+window->parentWidget()->baseSize().width();

      int actWidth = std::max(widthForEach, preferredWidth);

      window->parentWidget()->setGeometry( y, 0, actWidth , ws_->height() );
      y += actWidth;
    }
  }

  void TOPPViewBase::connectWindowSignals_(SpectrumWindow* sw)
  {
    connect(sw->widget()->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(updateToolbar()));
    connect(sw,SIGNAL(sendStatusMessage(std::string,OpenMS::UnsignedInt)),this,SLOT(showStatusMessage(std::string,OpenMS::UnsignedInt)));
    connect(sw,SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
    connect(sw,SIGNAL(modesChanged(QWidget*)),this,SLOT(updateToolbar()));
  
  	Spectrum2DWindow* sw2 = dynamic_cast<Spectrum2DWindow*>(sw);
  	if (sw2 != 0)
  	{
  		connect(sw2->getHorizontalProjection(),SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
  		connect(sw2->getVerticalProjection(),SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
  	}
  }

  //! returns selected peaks of active spectrum framed by \c layer_index_.begin() and the last peak BEFORE \c layer_index_.end();
  vector<MSExperiment<>::SpectrumType::Iterator> TOPPViewBase::getActiveSpectrumSelectedPeaks()
  {
    Spectrum1DWindow* w1 = active1DWindow_();
    if (w1)
    {
      return (w1->widget()->canvas()->getSelectedPeaks());
    }
    return vector<MSExperiment<>::SpectrumType::Iterator>();
  }

  void TOPPViewBase::gotoDialog()
  {
    SpectrumWindow* w = activeWindow_();
    if (w)
    {
      w->showGoToDialog();
    }
  }


  void TOPPViewBase::showPeaklistActiveSpectrum()
  {
    std::vector<MSExperiment<>::SpectrumType::Iterator> peaks = getActiveSpectrumSelectedPeaks();

    //only if peaks selected:
    if (peaks.size() >= 3)
    {
      std::ostringstream peaklist;

      for (std::vector<MSExperiment<>::SpectrumType::Iterator>::iterator it = peaks.begin(); it != peaks.end(); it++)
      {
        if (it == peaks.begin())
        {
          peaklist << "First Peak in Spectrum: Position " << (**it).getPosition()[0] << ", Intensity " << (**it).getIntensity() << endl;
        }
        else if (*it == *(peaks.rbegin()))
        {
          peaklist << "Last Peak in Spectrum: Position " << (**it).getPosition()[0] << ", Intensity " << (**it).getIntensity();
        }
        else
        {
          peaklist << "Selected Peak: Position " << (**it).getPosition()[0] << ", Intensity " << (**it).getIntensity() << endl;
        }
      }
      QMessageBox::information(this, QString("Peaklist"), QString(peaklist.str().c_str()));
    }
    else
    {
      QMessageBox::information(this, QString("Peaklist"), QString("No Peaks selected!"));
    }

  }

  void TOPPViewBase::smoothActiveSpectrum()
  {
    SmoothingDialog dialog(this);
    if (dialog.exec())
    {
      float kernel_width = dialog.getKernelWidth();
			float spacing = dialog.getSpacing();
			bool resampling_flag = dialog.getResampling();
      bool sgolay_flag = dialog.getSGolay();
      
      // 1D smoothing
      Spectrum1DWindow* w = active1DWindow_();
      if (w!=0)
      {
      	// raw data spectrum
        const Spectrum1DCanvas::ExperimentType& exp_raw = w->widget()->canvas()->getCurrentPeakData();
				
				//add new layer
        String new_name = w->widget()->canvas()->getCurrentLayer().name+" (smoothed)";
        Spectrum1DCanvas::ExperimentType& exp_smoothed = w->widget()->canvas()->addEmptyPeakLayer();

        // add one spectrum
        exp_smoothed.resize(1);
      	if (sgolay_flag)
        {
        	SavitzkyGolaySVDFilter sgolay;
        	try
	        {
	          sgolay.setOrder(dialog.getSGolayOrder());
	
	          // if resampling take the resampling spacing
	          if (resampling_flag)
	          { 
		          // and determine the number of kernel coefficients
		          int frame_size = (int) ceil(kernel_width / spacing + 1);
		
		          // the number has to be odd
		          if (!isOdd(frame_size))
		          {
		            frame_size += 1;
		          }
		          sgolay.setWindowSize(frame_size);
		          
		          Spectrum1DCanvas::ExperimentType::SpectrumType resampled_spectrum;
		          LinearResampler resampler;
		          resampler.setSpacing(spacing);
		          resampler.raster(exp_raw[0],resampled_spectrum);
		    
		          sgolay.filter(resampled_spectrum,exp_smoothed[0]);
		        }
		        // savitzky golay without resampling
		        else
		        {
		        	// else compute the spacing of data
			        float s =  ((exp_raw[0].end()-1)->getPos() - exp_raw[0].begin()->getPos()) / (exp_raw[0].size() + 1);
			        
			        // and determine the number of kernel coefficients
			        int frame_size = (int) ceil(kernel_width / s + 1);
			
			        // the number has to be odd
			        if (!isOdd(frame_size))
			        {
			         frame_size += 1;
			        }
			        sgolay.setWindowSize(frame_size);      
		          sgolay.filter(exp_raw[0],exp_smoothed[0]);
	          }
	        }
	        catch (Exception::InvalidValue)
	        {
	        	QMessageBox::warning(this,"Order of the smoothing polynomial too large!","Choose a smaller order or enlarge the width of the smoothing kernel.");
	        	w->widget()->canvas()->finishAdding();
	        	w->widget()->canvas()->setCurrentLayerName(new_name);
						return;
	        }
        }
        // gaussian smoothing
        else
        {
        	GaussFilter gauss;
        	
          gauss.setKernelWidth(kernel_width);
          gauss.filter(exp_raw[0],exp_smoothed[0]);
        }

        //color smoothed data
        for (Spectrum1DCanvas::ExperimentType::SpectrumType::Iterator it = exp_smoothed[0].begin(); it!= exp_smoothed[0].end(); ++it)
        {
          it->setMetaValue(UnsignedInt(5),string("#FF00FF"));
        }

        w->widget()->canvas()->finishAdding();
        w->widget()->canvas()->setCurrentLayerName(new_name);
      }
      // 2D window
      else
      {
        Spectrum2DWindow* w2 = active2DWindow_();
        if (w2!=0)
        {
          const Spectrum2DCanvas::ExperimentType& exp_raw = w2->widget()->canvas()->getCurrentPeakData();
          	
         	//add new window for picked peaks
          Spectrum2DWindow* w_smoothed = new Spectrum2DWindow(ws_);
					
          //set main preferences
          w_smoothed->setMainPreferences(prefs_);
          String new_name = w2->widget()->canvas()->getCurrentLayer().name+" (smoothed)";

          Spectrum2DCanvas::ExperimentType& exp2 = w_smoothed->widget()->canvas()->addEmptyPeakLayer();

          String filename = w2->widget()->canvas()->getCurrentLayer().name+"_smoothed";

					// savitzky golay filtering

					if (sgolay_flag)
        	{
	        	SavitzkyGolaySVDFilter sgolay;
	        	try
	        	{
		          sgolay.setOrder(dialog.getSGolayOrder());
	
		          // if resampling take the resampling spacing
		          if (resampling_flag)
		          { 
			          // and determine the number of kernel coefficients
			          int frame_size = (int) ceil(kernel_width / spacing + 1);
			
			          // the number has to be odd
			          if (!isOdd(frame_size))
			          {
			            frame_size += 1;
			          }
			          sgolay.setWindowSize(frame_size);
		          
			          Spectrum2DCanvas::ExperimentType resampled_experiment;
	            	LinearResampler lin_resampler;
	            	lin_resampler.setSpacing(spacing);
	
	            	unsigned int n = exp_raw.size();
	            	// resample and filter every scan
	            	for (unsigned int i = 0; i < n; ++i)
	            	{
		              // temporary container for the resampled data
		              Spectrum2DCanvas::ExperimentType::SpectrumType resampled_data;
		              lin_resampler.raster(exp_raw[i],resampled_data);
		
		              Spectrum2DCanvas::ExperimentType::SpectrumType spectrum;
		              sgolay.filter(resampled_data, spectrum);
	
		              // if any peaks are found copy the spectrum settings
		              if (spectrum.size() > 0)
		              {
	                	exp2.push_back(spectrum);
	              	}
	            	}
	          	}
	          	// savitzky golay filtering without resampling
	          	else
	          	{
	          		// else compute the spacing of data
				        float s =  ((exp_raw[0].end()-1)->getPos() - exp_raw[0].begin()->getPos()) / (exp_raw[0].size() + 1);
				        
				        // and determine the number of kernel coefficients
				        int frame_size = (int) ceil(kernel_width / s + 1);
				
				        // the number has to be odd
				        if (!isOdd(frame_size))
				        {
			         		frame_size += 1;
			        	}
				        sgolay.setWindowSize(frame_size);      
	          		sgolay.filterExperiment(exp_raw,exp2);
	            }
          	}
            catch (Exception::InvalidValue)
	        	{
	        		QMessageBox::warning(this,"Order of the smoothing polynomial too large!","Choose a smaller order or enlarge the width of the smoothing kernel.");
	        	
	        		delete w_smoothed;
							return;
	        	}
          }

          // gaussian filtering
          else
          {
           	GaussFilter gauss;
           	
           	gauss.setKernelWidth(kernel_width);
            gauss.filterExperiment(exp_raw,exp2);
          }
           
	      	w_smoothed->widget()->canvas()->finishAdding();
					w_smoothed->widget()->canvas()->setCurrentLayerName(new_name);
					
	        connectWindowSignals_(w_smoothed);
	        w_smoothed->setWindowTitle(new_name.c_str());
	        addClient(w_smoothed,filename);
	        addTab_(w_smoothed,new_name);
	
	        w_smoothed->showMaximized();
	        //           String gradient_peaks("Linear|0,#dbffcf;100,#00ff00");
	        //           w_smoothed->widget()->canvas()->setDotGradient(gradient_peaks);
	      }
	    }
      updateLayerbar();
    }
  }



  void TOPPViewBase::baselineFilteringActiveSpectrum()
  {
    BaselineFilteringDialog dialog(this);
    if (dialog.exec())
    {
      // 1D baseline filtering
      Spectrum1DWindow* w = active1DWindow_();
      if (w!=0)
      {
        const Spectrum1DCanvas::ExperimentType& exp_raw = w->widget()->canvas()->getCurrentPeakData();
        TopHatFilter tophat;
        tophat.setStrucElemSize(dialog.getStrucElemWidth());

        bool resampling_flag = dialog.getResampling();

        //add new layer
        String new_name = w->widget()->canvas()->getCurrentLayer().name+" (basline)";
        Spectrum1DCanvas::ExperimentType& exp_filtered = w->widget()->canvas()->addEmptyPeakLayer();

        // add one spectrum
        exp_filtered.resize(1);

        // Resampling before baseline filtering?
        if (resampling_flag)
        {
          Spectrum1DCanvas::ExperimentType::SpectrumType resampled_spectrum;
          LinearResampler resampler;
          resampler.setSpacing(dialog.getSpacing());
          resampler.raster(exp_raw[0],resampled_spectrum);
          tophat.filter(resampled_spectrum,exp_filtered[0]);
        }
        else
        {
          tophat.filter(exp_raw[0],exp_filtered[0]);
        }
        //color smoothed data
        for (Spectrum1DCanvas::ExperimentType::SpectrumType::Iterator it = exp_filtered[0].begin(); it!= exp_filtered[0].end(); ++it)
        {
          it->setMetaValue(UnsignedInt(5),string("#FF00FF"));
        }

        w->widget()->canvas()->finishAdding();
        w->widget()->canvas()->setCurrentLayerName(new_name);
      }
      else
      {
        Spectrum2DWindow* w2 = active2DWindow_();
        if (w2!=0)
        {
          const Spectrum2DCanvas::ExperimentType& exp_raw = w2->widget()->canvas()->getCurrentPeakData();
          TopHatFilter tophat;
          tophat.setStrucElemSize(dialog.getStrucElemWidth());

          bool resampling_flag = dialog.getResampling();

          //add new window for baseline filtered data
          Spectrum2DWindow* w_tophat = new Spectrum2DWindow(ws_);
          
          //set main preferences
          w_tophat->setMainPreferences(prefs_);
          String new_name = w2->widget()->canvas()->getCurrentLayer().name+" (basline)";

          Spectrum2DCanvas::ExperimentType& exp_filtered = w_tophat->widget()->canvas()->addEmptyPeakLayer();

          String filename = w2->widget()->canvas()->getCurrentLayer().name+"_filtered";

          // Resampling before baseline filtering?
          if (resampling_flag)
          {
            Spectrum2DCanvas::ExperimentType resampled_experiment;
            LinearResampler lin_resampler;
            lin_resampler.setSpacing(dialog.getSpacing());

            unsigned int n = exp_raw.size();
            // resample and filter every scan
            for (unsigned int i = 0; i < n; ++i)
            {
              // temporary container for the resampled data
              Spectrum2DCanvas::ExperimentType::SpectrumType resampled_data;
              lin_resampler.raster(exp_raw[i],resampled_data);

              Spectrum2DCanvas::ExperimentType::SpectrumType spectrum;
              tophat.filter(resampled_data, spectrum);

              // if any peaks are found copy the spectrum settings
              if (spectrum.size() > 0)
              {
                // copy the spectrum settings
                //  static_cast<SpectrumSettings&>(spectrum) = ms_exp_raw[i];
                // spectrum.setType(SpectrumSettings::RAWDATA);

                // copy the spectrum information
                //                 spectrum.getPrecursorPeak() = ms_exp_raw[i].getPrecursorPeak();
                //                 spectrum.setRetentionTime(ms_exp_raw[i].getRetentionTime());
                //                 spectrum.setMSLevel(ms_exp_raw[i].getMSLevel());
                //                 spectrum.getName() = ms_exp_raw[i].getName();

                exp_filtered.push_back(spectrum);
              }
            }
          }
          else
          {
              tophat.filterExperiment(exp_raw,exp_filtered);
          }

          w_tophat->widget()->canvas()->finishAdding();
          w_tophat->widget()->canvas()->setCurrentLayerName(new_name);
          connectWindowSignals_(w_tophat);
          w_tophat->setWindowTitle(new_name.c_str());
          addClient(w_tophat,filename);
          addTab_(w_tophat,new_name);

          w_tophat->showMaximized();

          //           String gradient_peaks("Linear|0,#dbffcf;100,#00ff00");
          //           w_smoothed->widget()->canvas()->setDotGradient(gradient_peaks);
        }
      }
      updateLayerbar();
    }
  }




  void TOPPViewBase::pickActiveSpectrum()
  {
    PeakPickingDialog dialog(this);
    if (dialog.exec())
    {
      //pick data
      PeakPickerCWT peak_picker;
      //estimate the scale a given the fwhm f
      //the fwhm of a marr wavelet is f ~ 1.252*a
      peak_picker.setWaveletScale(dialog.getFwhm() / 1.252);
      peak_picker.setFwhmBound(dialog.getFwhm());
      peak_picker.setPeakBound(dialog.getPeakHeight());
      peak_picker.setPeakBoundMs2Level(dialog.getPeakHeightMs2());
      peak_picker.setOptimizationFlag(dialog.getOptimization());
      peak_picker.setSignalToNoiseLevel(dialog.getSignalToNoise());

      Spectrum1DWindow* w = active1DWindow_();
      if (w!=0)
      {
        //add new layer
        String new_name = w->widget()->canvas()->getCurrentLayer().name+" (picked)";
        Spectrum1DCanvas::ExperimentType& exp = w->widget()->canvas()->addEmptyPeakLayer();

        // add one spectrum
        exp.resize(1);
        // pick peaks
        peak_picker.pick(w->widget()->canvas()->getCurrentPeakData()[0].begin(),w->widget()->canvas()->getCurrentPeakData()[0].end(),exp[0]);

        //color picked peaks
        for (Spectrum1DCanvas::ExperimentType::SpectrumType::Iterator it = exp[0].begin(); it!= exp[0].end(); ++it)
        {
          it->setMetaValue(UnsignedInt(5),string("#FF00FF"));
        }

        w->widget()->canvas()->finishAdding();
        w->widget()->canvas()->setCurrentLayerName(new_name);
      }
      else
      {

        Spectrum2DWindow* w2 = active2DWindow_();
        if (w2!=0)
        {
          //add new window for picked peaks
          Spectrum2DWindow* w_picked = new Spectrum2DWindow(ws_);
          //set main preferences
          w_picked->setMainPreferences(prefs_);
          String new_name = w2->widget()->canvas()->getCurrentLayer().name+"(picked)";

          Spectrum2DCanvas::ExperimentType& exp2 = w_picked->widget()->canvas()->addEmptyPeakLayer();

          String filename = w2->widget()->canvas()->getCurrentLayer().name+"_picked";

          peak_picker.pickExperiment(w2->widget()->canvas()->getCurrentPeakData(),exp2);

          w_picked->widget()->canvas()->finishAdding();
          w_picked->widget()->canvas()->setCurrentLayerName(new_name);

          connectWindowSignals_(w_picked);
          w_picked->setWindowTitle(new_name.c_str());
          addClient(w_picked,filename);
          addTab_(w_picked,new_name);

          w_picked->showMaximized();
        }
      }
      updateLayerbar();
    }
  }

  SpectrumWindow*  TOPPViewBase::activeWindow_() const
  {
    return dynamic_cast<SpectrumWindow*>(ws_->activeWindow());
  }

  SpectrumCanvas*  TOPPViewBase::activeCanvas_() const
  {
    SpectrumWindow* sw = dynamic_cast<SpectrumWindow*>(ws_->activeWindow());
    if (sw == 0)
    {
    	return 0;
    }
    return sw->widget()->canvas();
  }

  Spectrum1DWindow* TOPPViewBase::active1DWindow_() const
  {
    Spectrum1DWindow* s1;
    if ((s1 = dynamic_cast<Spectrum1DWindow*>(ws_->activeWindow())))
    {
      return s1;
    }
    return 0;
  }

  Spectrum2DWindow* TOPPViewBase::active2DWindow_() const
  {
    Spectrum2DWindow* s2;
    if ((s2 = dynamic_cast<Spectrum2DWindow*>(ws_->activeWindow())))
    {
      return s2;
    }
    return 0;
  }

  Spectrum3DWindow* TOPPViewBase::active3DWindow_() const
  {
    Spectrum3DWindow* s3;
    if ((s3 = dynamic_cast<Spectrum3DWindow*>(ws_->activeWindow())))
    {
      return s3;
    }
    return 0;
  }

  PreferencesDialogPage* TOPPViewBase::createPreferences(QWidget* parent)
  {
    return new TOPPViewBasePDP(this,parent);
  }

  void TOPPViewBase::loadPreferences(string filename)
  {
    //compose default ini file path
    String default_ini_file;
    char * home;
    home = getenv ("HOME");
    if (home!=NULL)
    {
      default_ini_file = home;
      default_ini_file = default_ini_file + "/";
    }
    default_ini_file = default_ini_file + ".TOPPView.ini";

    if (filename=="")
    {
      filename = default_ini_file;
    }
    prefs_.setValue("PreferencesFile" , filename);

    //load preferences, if file exists
    if (File::exists(filename))
    {
      prefs_.clear();
      prefs_.setValue("PreferencesFile" , filename);
      prefs_.load(filename);
    }
    else
    {
      if (filename != default_ini_file)
      {
        cerr << "Unable to load INI File: '" << filename << "'" << endl;
      }
    }

    //set missing defaults
    checkPreferences_();

    //set the recent files
    Param p = prefs_.copy("Preferences:RecentFiles");
    if (p.size()!=0)
    {
      for (Param::ConstIterator it=p.begin() ; it!=p.end() ; ++it)
      {
      	QString filename = ((string)it->second).c_str();
      	if (File::exists(filename.toAscii().data())) recent_files_.append(filename);
      }
    }
		
    updateRecentMenu_();
  }

  void TOPPViewBase::savePreferences()
  {
    // replace recent files
    prefs_.remove("Preferences:RecentFiles");

    for (int i = 0; i < recent_files_.size(); ++i)
    {
      prefs_.setValue("Preferences:RecentFiles:"+String(i),recent_files_[i].toStdString());
    }

    //save only the subsection that begins with "Preferences:"
    try
    {
      prefs_.copy("Preferences:").store(string(prefs_.getValue("PreferencesFile")));
    }
    catch(Exception::UnableToCreateFile& e)
    {
      cerr << "Unable to create INI File: '" << string(prefs_.getValue("PreferencesFile")) << "'" << endl;
    }
  }

  void TOPPViewBase::checkPreferences_()
  {
    Param default_preferences;

    //general
    default_preferences.setValue("DB:Host", "localhost");
    default_preferences.setValue("DB:Login", "NoName");
    default_preferences.setValue("DB:Name", "OpenMS");
    default_preferences.setValue("DB:Port", "3306");
    default_preferences.setValue("DefaultMapView", "2D");
    default_preferences.setValue("DefaultPath", ".");
    default_preferences.setValue("NumberOfRecentFiles", 15);
    default_preferences.setValue("Legend", "Hide");
    default_preferences.setValue("MapIntensityCutoff", "None");

    //1d
    default_preferences.setValue("1D:HighColor", "#ff0000");
    default_preferences.setValue("1D:IconColor", "#000000");
    default_preferences.setValue("1D:PeakColor", "#0000ff");
    default_preferences.setValue("1D:BackgroundColor", "#ffffff");
    default_preferences.setValue("1D:Mapping:MappingOfMzTo","X-Axis");

    //2d
    default_preferences.setValue("2D:BackgroundColor", "#ffffff");
    default_preferences.setValue("2D:MarchingSquaresSteps", 20);
    default_preferences.setValue("2D:InterpolationSteps", 200);
    default_preferences.setValue("2D:Dot:Gradient", "Linear|0,#efef00;11,#ffaa00;32,#ff0000;55,#aa00ff;78,#5500ff;100,#000000");
    default_preferences.setValue("2D:Dot:Mode", 1);
    default_preferences.setValue("2D:Surface:Gradient", "Linear|0,#ffffff;22,#fdffcb;50,#ffb4b4;75,#d7cfff;100,#c1c1c1");
    default_preferences.setValue("2D:Contour:Lines", 8);
    default_preferences.setValue("2D:Mapping:MappingOfMzTo","X-Axis");

    //3d
    default_preferences.setValue("3D:Dot:Mode", 1);
    default_preferences.setValue("3D:Shade:Mode", 1);
    default_preferences.setValue("3D:Dot:Gradient", "Linear|0,#efef00;11,#ffaa00;32,#ff0000;55,#aa00ff;78,#5500ff;100,#000000");
    default_preferences.setValue("3D:Dot:InterpolationSteps",200);
    default_preferences.setValue("3D:BackgroundColor", "#ffffff");
    default_preferences.setValue("3D:AxesColor", "#000000");
    default_preferences.setValue("3D:IntMode",0);
    default_preferences.setValue("3D:Dot:LineWidth",2);
    default_preferences.setValue("3D:IntScale:Mode",0);
		default_preferences.setValue("3D:Data:Mode",0);
		default_preferences.setValue("3D:Data:Reduction:Max",10);
		default_preferences.setValue("3D:Data:Reduction:Sum",10);
		default_preferences.setValue("3D:DisplayedPeaks",10000);
		default_preferences.setValue("3D:Reduction:Mode","MaxReduction");

    prefs_.setDefaults(default_preferences,"Preferences");
  }

  void TOPPViewBase::openRecentFile()
  {
  	//cout << "TEST" << endl;

		QAction* action = qobject_cast<QAction *>(sender());
    if (action)
		{
	  	String filename = action->data().toString().toStdString();
	  	//cout << "filename: " << filename << endl;
	  	setCursor(Qt::WaitCursor);
 	  	OpenDialog::Mower mow = OpenDialog::NO_MOWER;
			if ( getPrefAsString("Preferences:MapIntensityCutoff")=="Noise Estimator")
			{
				mow = OpenDialog::NOISE_ESTIMATOR;
			}
   		addSpectrum(filename,true,getPrefAsString("Preferences:DefaultMapView")=="2D",true,mow);
			setCursor(Qt::ArrowCursor); 	
			//cout << "TEST_END" << endl;
		}
 }

  void TOPPViewBase::findFeaturesActiveSpectrum()
  {
#ifdef CGAL_DEF
    Spectrum2DWindow* w = active2DWindow_();
    if (w!=0)
    {
      FeaFiDialog dialog(this);
      if (dialog.exec() == QDialog::Accepted)
      {
        //find features
        FeatureFinder& finder = dialog.getFeatureFinder();
        
        finder.setData(w->widget()->canvas()->getCurrentPeakData().begin(),w->widget()->canvas()->getCurrentPeakData().end(),1500);
        DFeatureMap<2> map = finder.run();

        //display features
        w->widget()->canvas()->addLayer(map,false);
        w->widget()->canvas()->setCurrentLayerName(w->widget()->canvas()->getCurrentLayer().name+" (features)");
        updateLayerbar();
      }
    }
#endif
  }

  void TOPPViewBase::openSpectrumDialog()
  {
    OpenDialog dialog(prefs_,this);
    if (dialog.exec())
    {
      //Open Files
      setCursor(Qt::WaitCursor);
      if (dialog.getSource()==OpenDialog::FILE)
      {
        for(vector<String>::const_iterator it=dialog.getNames().begin();it!=dialog.getNames().end();it++)
        {
          addSpectrum(*it,dialog.isOpenAsNewTab(),dialog.isViewMaps2D(),true,dialog.getMower(),dialog.forcedFileType());
        }
      }
      else
        // Open from DB
      {
        for(vector<String>::const_iterator it=dialog.getNames().begin();it!=dialog.getNames().end();it++)
        {
          addDBSpectrum(it->toInt(),dialog.isOpenAsNewTab(),dialog.isViewMaps2D(),true,dialog.getMower());
        }
      }
      setCursor(Qt::ArrowCursor);
      maximizeActiveSpectrum();
    }
  }

} //namespace OpenMS


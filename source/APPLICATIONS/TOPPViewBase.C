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
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DB/DBAdapter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/VISUAL/DIALOGS/SaveImageDialog.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum3DWidget.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureMapFile.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolaySVDFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/BASELINE/TopHatFilter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>
#include <OpenMS/FORMAT/FeaturePairsFile.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewPrefDialog.h>
#include <OpenMS/ANALYSIS/ID/IDFeatureMapper.h>
#include <OpenMS/ANALYSIS/ID/IDSpectrumMapper.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

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
#include <QtCore/QDate>
#include <QtGui/QWhatsThis>

//action modes
#include "../VISUAL/ICONS/zoom.xpm"
#include "../VISUAL/ICONS/selection.xpm"

//intensity modes
#include "../VISUAL/ICONS/lin.xpm"
#include "../VISUAL/ICONS/logarithm.xpm"
#include "../VISUAL/ICONS/percentage.xpm"
#include "../VISUAL/ICONS/snap.xpm"

//common
#include "../VISUAL/ICONS/reset_zoom.xpm"
#include "../VISUAL/ICONS/tile_horizontal.xpm"
#include "../VISUAL/ICONS/tile_vertical.xpm"

//1d
#include "../VISUAL/ICONS/lines.xpm"
#include "../VISUAL/ICONS/peaks.xpm"

//2d
#include "../VISUAL/ICONS/colors.xpm"
#include "../VISUAL/ICONS/contours.xpm"
#include "../VISUAL/ICONS/precursors.xpm"
#include "../VISUAL/ICONS/projections.xpm"
#include "../VISUAL/ICONS/convexhulls.xpm"
#include "../VISUAL/ICONS/numbers.xpm"

#include "../VISUAL/ICONS/TOPPView.xpm"

#include <algorithm>

using namespace std;

namespace OpenMS
{
  using namespace Internal;
	using namespace Math;

  TOPPViewBase::TOPPViewBase(QWidget* parent):
      QMainWindow(parent),
      DefaultParamHandler("TOPPViewBase")
  {
  	setWindowTitle("TOPPView");
    setWindowIcon(QIcon(toppview));
    //prevents errors caused by too small width,height values
    setMinimumSize(400,400);

    // create dummy widget (to be able to have a layout), Tab bar and workspace
    QWidget* dummy = new QWidget(this);
    setCentralWidget(dummy);
    QVBoxLayout* box_layout = new QVBoxLayout(dummy);
    tab_bar_ = new EnhancedTabBar(dummy);
    tab_bar_->setWhatsThis("Double-click tab to close it.");
    tab_bar_->addTab("dummy");
    tab_bar_->setMinimumSize(tab_bar_->sizeHint());
    tab_bar_->removeTab(0);

    //connect slots and sigals for selecting spectra
    connect(tab_bar_,SIGNAL(currentChanged(int)),this,SLOT(focusByTab(int)));
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
    file->addAction("&Open",this,SLOT(openSpectrumDialog()), Qt::CTRL+Qt::Key_O);
    file->addAction("&Close",this,SLOT(closeFile()));
    file->addSeparator();
    file->addAction("&Edit INI file", this, SLOT(editParamDialog()));
		file->addSeparator();
    
    QMenu* recent_menu = new QMenu("&Recent files", this);
    recent_as_new_layer_ = recent_menu->addAction("open as new layer");
    recent_as_new_layer_->setCheckable(true);
    recent_menu->addSeparator();
    //create the max mumber of recent files actions
  	recent_actions_.resize(20);
		for (UInt i = 0; i<20; ++i)
		{
			recent_actions_[i] = recent_menu->addAction("",this,SLOT(openRecentFile()));
			recent_actions_[i]->setVisible(false);
		}
  	file->addMenu(recent_menu);

    file->addSeparator();
    file->addAction("&Preferences",this, SLOT(preferencesDialog()));
    file->addAction("&Quit",qApp,SLOT(quit()), Qt::CTRL+Qt::Key_Q);
    
    //Layer menu
    QMenu* layer = new QMenu("&Layer",this);
    menuBar()->addMenu(layer);
    layer->addAction("&Save visible data",this,SLOT(saveLayer()), Qt::CTRL+Qt::Key_S);
    layer->addAction("&Edit metadata",this,SLOT(editMetadata()));
    layer->addAction("&Intensity distribution",this,SLOT(layerIntensityDistribution()));
		layer->addSeparator();
    layer->addAction("Apply &TOPP tool", this, SLOT(showTOPPDialog()), Qt::CTRL+Qt::Key_T);
    layer->addAction("&Annotate with identifiction", this, SLOT(annotateWithID()), Qt::CTRL+Qt::Key_A);
		layer->addSeparator();
    layer->addAction("&Preferences",this, SLOT(layerPreferencesDialog()));
    
    //View menu
    QMenu * view = new QMenu("&View",this);
    menuBar()->addMenu(view);
    view->addAction("&Go to",this,SLOT(gotoDialog()), Qt::CTRL+Qt::Key_G);
   	view->addAction("Show/Hide &axis legends",this,SLOT(changeAxisVisibility()));
   	view->addAction("Show/Hide &grid lines",this,SLOT(changeGridLines()));
   	
    //Image menu
    QMenu * image = new QMenu("&Image",this);
    menuBar()->addMenu(image);
    image->addAction("&Save to file",this,SLOT(saveImage()));
    image->addAction("&Print",this,SLOT(print()), Qt::CTRL+Qt::Key_P);

    //Windows menu
    QMenu * windows = new QMenu("&Windows", this);
    menuBar()->addMenu(windows);
    windows->addAction("&Cascade",this->ws_,SLOT(cascade()));
    windows->addAction("&Tile automatic",this->ws_,SLOT(tile()));
    windows->addAction(QIcon(QPixmap(tile_h)),"Tile &vertical",this,SLOT(tileHorizontal()));
    windows->addAction(QIcon(QPixmap(tile_v)),"Tile &horizontal",this,SLOT(tileVertical()));

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

    //set defaults
    //general
    defaults_.setValue("Preferences:DefaultMapView", "2D");
    defaults_.setValue("Preferences:DefaultPath", ".");
    defaults_.setValue("Preferences:TmpPath", "/tmp/");
    defaults_.setValue("Preferences:NumberOfRecentFiles", 15);
    defaults_.setValue("Preferences:Legend", "Show");
    defaults_.setValue("Preferences:MapIntensityCutoff", "None");
    //db
    defaults_.setValue("Preferences:DB:Host", "localhost");
    defaults_.setValue("Preferences:DB:Login", "NoName");
    defaults_.setValue("Preferences:DB:Name", "OpenMS");
    defaults_.setValue("Preferences:DB:Port", "3306");
    //1d
    defaults_.setValue("Preferences:1D:HighColor", "#ff0000");
    defaults_.setValue("Preferences:1D:IconColor", "#000000");
    defaults_.setValue("Preferences:1D:PeakColor", "#0000ff");
    defaults_.setValue("Preferences:1D:BackgroundColor", "#ffffff");
    //2d
    defaults_.setValue("Preferences:2D:BackgroundColor", "#ffffff");
    defaults_.setValue("Preferences:2D:MarchingSquaresSteps", 20);
    defaults_.setValue("Preferences:2D:InterpolationSteps", 200);
    defaults_.setValue("Preferences:2D:Dot:Gradient", "Linear|0,#efef00;7,#ffaa00;15,#ff0000;27,#aa00ff;55,#5500ff;100,#000000");
    defaults_.setValue("Preferences:2D:Surface:Gradient", "Linear|0,#ffffff;7,#fdffcb;20,#ffb4b4;50,#d7cfff;100,#c1c1c1");
    defaults_.setValue("Preferences:2D:Contour:Lines", 8);
    defaults_.setValue("Preferences:2D:Mapping:MappingOfMzTo","X-Axis");
    //3d
    defaults_.setValue("Preferences:3D:Dot:ShadeMode", 1);
    defaults_.setValue("Preferences:3D:Dot:Gradient", "Linear|0,#efef00;11,#ffaa00;32,#ff0000;55,#aa00ff;78,#5500ff;100,#000000");
    defaults_.setValue("Preferences:3D:Dot:InterpolationSteps",200);
    defaults_.setValue("Preferences:3D:BackgroundColor", "#ffffff");
    defaults_.setValue("Preferences:3D:AxesColor", "#000000");
    defaults_.setValue("Preferences:3D:Dot:LineWidth",2);
		defaults_.setValue("Preferences:3D:DisplayedPeaks",10000);
		defaults_.setValue("Preferences:3D:ReductionMode","Max reduction");
		
		defaults_.setValue("Preferences:Version","none");
		
		subsections_.push_back("Preferences:RecentFiles");
		
  	defaultsToParam_();
  	
  	//load param file
    loadPreferences();    
  }

  TOPPViewBase::~TOPPViewBase()
  {
  	//cout << "DEST TOPPViewBase" << endl;
  	savePreferences();
  }

  void TOPPViewBase::closeEvent(QCloseEvent* /*event*/)
  {
  	ws_->closeAllWindows();
  }

  void TOPPViewBase::addDBSpectrum(UInt db_id, bool as_new_window, bool maps_as_2d, bool maximize, OpenDialog::Mower use_mower)
  {
    //DBConnection for all DB queries
    DBConnection con;
    con.connect(param_.getValue("Preferences:DB:Name"), param_.getValue("Preferences:DB:Login"),param_.getValue("DBPassword"),param_.getValue("Preferences:DB:Host"),(Int)param_.getValue("Preferences:DB:Port"));

    //DB adapter
    DBAdapter dba(con);

    String db_id_string(db_id);
    QSqlQuery result;
    con.executeQuery("SELECT count(id) from DATA_Spectrum where fid_MSExperiment='"+db_id_string+"' and MSLevel='1'",result);

    //tab caption
    String caption = "DB ("+db_id_string+")";

    SpectrumWidget* w;

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
        w = new Spectrum1DWidget(param_, ws_);

        //determine Spectrum id
        con.executeQuery("SELECT id from DATA_Spectrum where fid_MSExperiment='"+db_id_string+"' and MSLevel='1'",result);
        UID spectrum_id = result.value(0).toInt();

        //load data
        exp = &(w->canvas()->addEmptyPeakLayer());
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
          w = new Spectrum2DWidget(param_, ws_);

          //load spectrum
          exp = &(w->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, (*exp));
        }
        //create 3D view
        else
        {
          //cout << "NEW 3D" << endl;
        	// create 3D window
          w = new Spectrum3DWidget(param_, ws_);

          //load data
          exp = &(w->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, *exp);
        }
      }
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
          QMessageBox::critical(this,"Wrong file type",("You cannot open 1D data ("+db_id_string+") in a 2D window!<BR>Please open the file in new tab.").c_str());
          return;
        }
        else //open it
        {
          //determine Spectrum id
          con.executeQuery("SELECT id from DATA_Spectrum where fid_MSExperiment='"+db_id_string+"' and MSLevel='1'", result);
          UID spectrum_id = result.value(0).toInt();

          //load data
	        exp = &(w->canvas()->addEmptyPeakLayer());
	        exp->resize(1);
	        dba.loadSpectrum(spectrum_id,(*exp)[0]);
        }
      }
      //create 2D/3D view
      else
      {
        Spectrum2DWidget* w2 = active2DWindow_();
        Spectrum3DWidget* w3 = active3DWindow_();
        //wrong active window type
        if (w2==0 && w3==0)
        {
          QMessageBox::critical(this,"Wrong file type",("You cannot open 1D data ("+db_id_string+") in a 2D/3D window!<BR>Please open the file in new tab.").c_str());
          return;
        }
        //create 2D view
        if (w2!=0)
        {
        	//cout << "ACTIVE 2D" << endl;
          w = w2;

          //load data
          exp = &(w->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, (*exp));
        }
        //create 3D view
        else
        {
        	//cout << "ACTIVE 3D" << endl;
          w = w3;

          //load data
          exp = &(w->canvas()->addEmptyPeakLayer());
          dba.loadExperiment(db_id, (*exp));
        }
      }
		}
		
    //noise estimator
    float cutoff = 0;
    if(use_mower!=OpenDialog::NO_MOWER && exp->size()>1)
    {
      cutoff = estimateNoise_(*exp);
    }
    w->canvas()->finishAdding(cutoff);
		w->canvas()->setCurrentLayerName(caption);
    //use_mower

    //do for all windows
    if (as_new_window)
    {
      showAsWindow_(w,caption);
    }

    //do for all (in active and in new window, 1D/2D/3D)
    if(maximize)
    {
      w->showMaximized();
    }

    //do for all windows
    updateLayerbar();
  }

  float TOPPViewBase::estimateNoise_(const SpectrumCanvas::ExperimentType& exp)
  {
  	//test if no scans with MS-level 1 exist => prevent deadlock
  	bool ms1_present = false;
  	for (UInt i = 0; i < exp.size(); ++i)
  	{
  		if (exp[i].getMSLevel()==1)
  		{
  			ms1_present = true;
  			break;
  		}
  	}
  	if (!ms1_present)
  	{
  		return 0.0;
  	}
  	
    float noise = 0.0;
    UInt count = 0;
    srand(time(0));
    //cout << "size: " << exp.size() << endl;
    while (count<10)
    {
      UInt scan = (UInt)( (double)rand() / ((double)(RAND_MAX)+1.0f) * (double)(exp.size()-1) );

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
        //cout << "scan: "<< scan <<" Groesse: " << tmp.size() << " Index: " << (UInt)ceil((float)(tmp.size()-1)/1.25f) << " Wert: "<< tmp[(UInt)ceil((float)(tmp.size()-1)/1.25f)] << endl;
        noise += tmp[(UInt)ceil((float)(tmp.size()-1)/1.25f)];

        ++count;
      }
    }
    noise = noise / 10.0f;
    return noise;
  }

  void TOPPViewBase::preferencesDialog()
  {
		Internal::TOPPViewPrefDialog dlg(this);
		
		//get pointers
		QLineEdit* default_path = dlg.findChild<QLineEdit*>("default_path");
		QLineEdit* temp_path = dlg.findChild<QLineEdit*>("temp_path");
		QSpinBox* recent_files = dlg.findChild<QSpinBox*>("recent_files");
		QComboBox* map_default = dlg.findChild<QComboBox*>("map_default");
		QComboBox* map_cutoff = dlg.findChild<QComboBox*>("map_cutoff");
		
		QLineEdit* db_host = dlg.findChild<QLineEdit*>("db_host");
		QSpinBox* db_port = dlg.findChild<QSpinBox*>("db_port");
		QLineEdit* db_name = dlg.findChild<QLineEdit*>("db_name");
		QLineEdit* db_login = dlg.findChild<QLineEdit*>("db_login");
		
		ColorSelector* color_1D = dlg.findChild<ColorSelector*>("color_1D");
		ColorSelector* selected_1D = dlg.findChild<ColorSelector*>("selected_1D");
		ColorSelector* icon_1D = dlg.findChild<ColorSelector*>("icon_1D");

		MultiGradientSelector* peak_2D = dlg.findChild<MultiGradientSelector*>("peak_2D");
		MultiGradientSelector* surface_2D = dlg.findChild<MultiGradientSelector*>("surface_2D");
		QSpinBox* contours_2D  = dlg.findChild<QSpinBox*>("contours_2D");
		QSpinBox* marching_cells_2D  = dlg.findChild<QSpinBox*>("marching_cells_2D");
		QComboBox* mapping_2D = dlg.findChild<QComboBox*>("mapping_2D");

		MultiGradientSelector* peak_3D = dlg.findChild<MultiGradientSelector*>("peak_3D");
		QComboBox* shade_3D = dlg.findChild<QComboBox*>("shade_3D");
		QSpinBox* line_width_3D  = dlg.findChild<QSpinBox*>("line_width_3D");
		QComboBox* reduction_3D = dlg.findChild<QComboBox*>("reduction_3D");
		QSpinBox* reduction_peaks_3D  = dlg.findChild<QSpinBox*>("reduction_peaks_3D");
		
		//set General values
		default_path->setText(param_.getValue("Preferences:DefaultPath").toQString());
		temp_path->setText(param_.getValue("Preferences:TmpPath").toQString());
		recent_files->setValue((Int)param_.getValue("Preferences:NumberOfRecentFiles"));
		map_default->setCurrentIndex(map_default->findText(param_.getValue("Preferences:DefaultMapView").toQString()));
		map_cutoff->setCurrentIndex(map_cutoff->findText(param_.getValue("Preferences:MapIntensityCutoff").toQString()));		

		db_host->setText(param_.getValue("Preferences:DB:Host").toQString());
		db_port->setValue((Int)param_.getValue("Preferences:DB:Port"));
		db_name->setText(param_.getValue("Preferences:DB:Name").toQString());
		db_login->setText(param_.getValue("Preferences:DB:Login").toQString());
		
		color_1D->setColor(QColor(param_.getValue("Preferences:1D:PeakColor").toQString()));
		selected_1D->setColor(QColor(param_.getValue("Preferences:1D:HighColor").toQString()));
		icon_1D->setColor(QColor(param_.getValue("Preferences:1D:IconColor").toQString()));

		peak_2D->gradient().fromString(param_.getValue("Preferences:2D:Dot:Gradient"));
		surface_2D->gradient().fromString(param_.getValue("Preferences:2D:Surface:Gradient"));
		contours_2D->setValue(UInt(param_.getValue("Preferences:2D:Contour:Lines")));
		marching_cells_2D->setValue(UInt(param_.getValue("Preferences:2D:MarchingSquaresSteps")));
		mapping_2D->setCurrentIndex(mapping_2D->findText(param_.getValue("Preferences:2D:Mapping:MappingOfMzTo").toQString()));

		peak_3D->gradient().fromString(param_.getValue("Preferences:3D:Dot:Gradient"));
		shade_3D->setCurrentIndex((Int)param_.getValue("Preferences:3D:Dot:ShadeMode"));
		line_width_3D->setValue((Int)param_.getValue("Preferences:3D:Dot:LineWidth"));
		reduction_3D->setCurrentIndex(reduction_3D->findText(param_.getValue("Preferences:3D:ReductionMode").toQString()));	
		reduction_peaks_3D->setValue((Int)param_.getValue("Preferences:3D:DisplayedPeaks"));
		
		//execute dialog
		if (dlg.exec())
		{
			//load data to param
			param_.setValue("Preferences:DefaultPath", default_path->text().toAscii().data());
			param_.setValue("Preferences:TmpPath", temp_path->text().toAscii().data());
			param_.setValue("Preferences:NumberOfRecentFiles", recent_files->value());
			param_.setValue("Preferences:DefaultMapView", map_default->currentText().toAscii().data());
			param_.setValue("Preferences:MapIntensityCutoff", map_cutoff->currentText().toAscii().data());

			param_.setValue("Preferences:DB:Host",db_host->text().toAscii().data());
			param_.setValue("Preferences:DB:Port",db_port->value());
			param_.setValue("Preferences:DB:Name",db_name->text().toAscii().data());
			param_.setValue("Preferences:DB:Login",db_login->text().toAscii().data());
			param_.remove("DBPassword");

			param_.setValue("Preferences:1D:PeakColor",color_1D->getColor().name().toAscii().data());
			param_.setValue("Preferences:1D:HighColor",selected_1D->getColor().name().toAscii().data());
			param_.setValue("Preferences:1D:IconColor",icon_1D->getColor().name().toAscii().data());

			param_.setValue("Preferences:2D:Dot:Gradient",peak_2D->gradient().toString());
			param_.setValue("Preferences:2D:Surface:Gradient",surface_2D->gradient().toString());
			param_.setValue("Preferences:2D:Contour:Lines",contours_2D->value());
			param_.setValue("Preferences:2D:MarchingSquaresSteps",marching_cells_2D->value());
			param_.setValue("Preferences:2D:Mapping:MappingOfMzTo",mapping_2D->currentText().toAscii().data());

			param_.setValue("Preferences:3D:Dot:Gradient",peak_3D->gradient().toString());
			param_.setValue("Preferences:3D:Dot:ShadeMode", shade_3D->currentIndex());
			param_.setValue("Preferences:3D:Dot:LineWidth",line_width_3D->value());
			param_.setValue("Preferences:3D:ReductionMode", reduction_3D->currentText().toAscii().data());
			param_.setValue("Preferences:3D:DisplayedPeaks",	reduction_peaks_3D->value());

			savePreferences();
		}
  }

  void TOPPViewBase::addSpectrum(const String& filename,bool as_new_window, bool maps_as_2d, bool maximize, OpenDialog::Mower use_mower, FileHandler::Type force_type)
  {
    if (!File::exists(filename))
    {
      QMessageBox::critical(this,"Open file error",("The file '"+filename+"' does not exist!").c_str());
      return;
    }

    //extract the filename without path
    String caption = File::basename(filename);

    //update recent files list
    addRecentFile_(filename);

    //windowpointer
    SpectrumWidget* w=0;
    
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
      QMessageBox::critical(this,"Open file error",("Could not determine file type of '"+filename+"'!").c_str());
      return;
		}

#ifdef DEBUG_TOPP
		cout << "TOPPViewBase::addSpectrum():";
		cout << " - File name: " << filename << endl;
		cout << " - File type: " << fh.typeToName(force_type) << endl;
		cout << " - New Window: " << as_new_window << endl;
		cout << " - Map as 2D: " << maps_as_2d << endl;		
#endif

    if (as_new_window)
    {
      if (force_type==FileHandler::DTA)
      {
        w = new Spectrum1DWidget(param_, ws_);
      }
      else if (maps_as_2d || force_type==FileHandler::FEATURE || force_type==FileHandler::FEATURE_PAIRS) //2d or features
      {
        w = new Spectrum2DWidget(param_, ws_);
      }
      else //3d
      {
        w = new Spectrum3DWidget(param_, ws_);
      }
    }
    else //!as_new_window
    {
      if (active1DWindow_()!=0) //1D window
      {
        w = active1DWindow_();
      }
      else if (active2DWindow_()!=0) //2d window
      {
        if (force_type==FileHandler::DTA)
        {
          QMessageBox::critical(this,"Wrong file type",("You cannot open 1D data ("+filename+") in a 2D window!").c_str());
          return;
        }
        w = active2DWindow_();
      }
      else if (active3DWindow_()!=0)//3d window
      {
        w = active3DWindow_();
      }
    }

    //try to read the data from file
    if (force_type==FileHandler::FEATURE) //features
    {
      FeatureMap<> map;
      try
      {
        FeatureMapFile().load(filename,map);
      }
      catch(Exception::Base& e)
      {
        QMessageBox::critical(this,"Error",(String("Error while reading feature file: ")+e.what()).c_str());
        return;
      }
      w->canvas()->addLayer(map,false);
      w->canvas()->setCurrentLayerName(caption);
    }
    else if (force_type==FileHandler::FEATURE_PAIRS) //feature pairs
    {
    	//load pairs
      std::vector< ElementPair < Feature > >  pairs;
      try
      {
        FeaturePairsFile().load(filename,pairs);
      }
      catch(Exception::Base& e)
      {
        QMessageBox::critical(this,"Error",(String("Error while reading feature pairs file: ")+e.what()).c_str());
        return;
      }
      
      //convert to features
      FeatureMap<> map;
      FeaturePairsFile::pairsToFeatures(pairs,map);
      w->canvas()->addLayer(map,true);
      w->canvas()->setCurrentLayerName(caption);
    }
    else
    {
      //try to read the data from file (raw/peak data)
      SpectrumCanvas::ExperimentType* exp = &(w->canvas()->addEmptyPeakLayer());
      try
      {
        fh.loadExperiment(filename,*exp, force_type);
      }
      catch(Exception::Base& e)
      {
        QMessageBox::critical(this,"Error",(String("Error while reading data file: ")+e.what()).c_str());
        return;
      }

      //check if only one scan is in a 2d file
      if (as_new_window && active1DWindow_()==0 && exp->size()==1)
      {
        delete(w);
        w = new Spectrum1DWidget(param_, ws_);
        exp = &(w->canvas()->addEmptyPeakLayer());
        FileHandler().loadExperiment(filename,*exp, force_type);
      }

      //do for all (in active and in new window, 1D/2D/3D)
      float cutoff=0;

      //calculate noise
      if(use_mower!=OpenDialog::NO_MOWER && exp->size()>1)
      {
        cutoff = estimateNoise_(*exp);
      }
      w->canvas()->finishAdding(cutoff);
      w->canvas()->setCurrentLayerName(caption);
    }
    
  	updateLayerbar();

    if (as_new_window)
    {
      showAsWindow_(w,caption);
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
		UInt number_of_recent_files = UInt(param_.getValue("Preferences:NumberOfRecentFiles"));
		while ((UInt)recent_files_.size() > number_of_recent_files)
		{
			recent_files_.removeLast();
		}
		
    updateRecentMenu_();
  }

  void TOPPViewBase::updateRecentMenu_()
  {
    //get/correct number of recent files
		UInt number_of_recent_files = UInt(param_.getValue("Preferences:NumberOfRecentFiles"));
		if (number_of_recent_files>20)
		{
			number_of_recent_files = 20;
			param_.setValue("Preferences:NumberOfRecentFiles",20);
		}
		
		for (UInt i = 0; i < 20; ++i)
		{
			if (i < (UInt)(recent_files_.size()))
			{
				recent_actions_[i]->setText(recent_files_[i]);
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

  void TOPPViewBase::addTab_(SpectrumWidget* w, const String& tabCaption)
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

  SpectrumWidget* TOPPViewBase::window_(int id) const
  {
  	//cout << "Looking for tab with id: " << id << endl;
  	QList<QWidget*> windows = ws_->windowList();
		for(int i=0; i< windows.size(); ++i)
		{
			SpectrumWidget* window = dynamic_cast<SpectrumWidget*>(windows.at(i));
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
  	SpectrumWidget* window = window_(tab_bar_->tabData(index).toInt());
  	if (window)
  	{
  		window->close();
  	}
  }
 
  void TOPPViewBase::focusByTab(int index)
  {
  	SpectrumWidget* window = window_(tab_bar_->tabData(index).toInt());
  	if (window)
  	{
  		window->setFocus();
  	}
  }

  void TOPPViewBase::saveImage()
  {
    //check if there is a active window
    SpectrumWidget* window = activeWindow_();
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
					if (File::writable(file_name.toStdString()))
					{
          	QImage image = window->getImage(dialog->getXSize(),dialog->getYSize());
          	image.save(file_name,format.toAscii().data(),100);
        	}
        	else
        	{
        	  QMessageBox::critical(this,"Error writing file!",(String("Cannot write to '")+file_name.toStdString()	+"'!").c_str());
        	}
        }
      }
    }
  }

  void TOPPViewBase::print()
  {
#ifndef QT_NO_PRINTER
    SpectrumWidget* window = activeWindow_();
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
        QImage image = window->getImage(body.width(),body.height());
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
      const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
      //warn if hidden layer => wrong layer selected...
    	if (!layer.visible)
    	{
    		QMessageBox::warning(this,"Warning","The current layer is not visible!");
    	}
    	//Visible area
    	const SpectrumCanvas::AreaType& area = activeWindow_()->canvas()->getVisibleArea();
    	
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
				
				UInt i = 0;
    		for (LayerData::ExperimentType::ConstIterator it=begin; it!=end; ++it)
    		{
  				out[i].SpectrumSettings::operator=(*it);
  				out[i].setRT(it->getRT());
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
    		  QMessageBox::critical(this,"Error","The displayed region of the current layer is empty!");
    		  return;
    		}
    		//one scan => DTA
    		else if (out.size()==1)
    		{
		      QString file_name = QFileDialog::getSaveFileName(this, "Save file", param_.getValue("Preferences:DefaultPath").toString().c_str(),
					                    "DTA files (*.dta)");
					if (!file_name.isEmpty())
					{
					  DTAFile().store(file_name.toAscii().data(),out[0]);
					}
    		}
    		//more than one scan => MzData
    		else
    		{
		      QString file_name = QFileDialog::getSaveFileName(this, "Save file",  param_.getValue("Preferences:DefaultPath").toString().c_str(),
					                    "MzData files (*.mzData)");
					if (!file_name.isEmpty())
					{
					  MzDataFile().store(file_name.toAscii().data(),out);
					}
    		}
			}
			else //feature data
			{
    		//Extract selected visible data to out
    		LayerData::FeatureMapType out;
    		out.ExperimentalSettings::operator=(layer.features);
    		for (LayerData::FeatureMapType::ConstIterator it=layer.features.begin(); it!=layer.features.end(); ++it)
    		{
					if ( it->getIntensity() >= layer.min_int && 
							 it->getIntensity() <= layer.max_int &&
							 it->getRT() >= area.min()[1] &&
							 it->getRT() <= area.max()[1] &&
							 it->getMZ() >= area.min()[0] &&
							 it->getMZ() <= area.max()[0]
						 )
					{
						out.push_back(*it);
					}
  			}
    		//no extracted data
    		if (out.size()==0)
    		{
    		  QMessageBox::critical(this,"Error","The displayed region of the current layer is empty!");
    		  return;
    		}
    		else
    		{
		      QString file_name = QFileDialog::getSaveFileName(this, "Save file",  param_.getValue("Preferences:DefaultPath").toString().c_str(),
					                    "features files (*.feat)" );
					if (!file_name.isEmpty())
					{
					  FeatureMapFile().store(file_name.toAscii().data(),out);
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
      const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
      //warn if hidden layer => wrong layer selected...
    	if (!layer.visible)
    	{
    		QMessageBox::warning(this,"Warning","The current layer is not visible!");
    	}
			MSMetaDataExplorer dlg(true, this);
      dlg.setWindowTitle("Edit meta data");
			if (layer.type==LayerData::DT_PEAK) //peak data
    	{
    		dlg.visualize(const_cast<LayerData&>(layer).peaks);
				
    	}
    	else //feature data
    	{
    		dlg.visualize(const_cast<LayerData&>(layer).features);
    	}
      dlg.exec();
    }
  }

  void TOPPViewBase::layerPreferencesDialog()
  {
    if (ws_->activeWindow())
    {
			activeWindow_()->canvas()->showPreferencesDialog();
    }
  }

  void TOPPViewBase::layerIntensityDistribution()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      activeWindow_()->showIntensityDistribution();
    }  	
  }

  void TOPPViewBase::changeAxisVisibility()
  {
    //check if there is a active window
    if (ws_->activeWindow())
    {
      activeWindow_()->showLegend(!activeWindow_()->isLegendShown());
    }  	
  }

  void TOPPViewBase::createToolBars_()
  {
  	QToolButton* b;
  	
  	//--Basic tool bar for all views--
    tool_bar_ = addToolBar("Basic tool bar");
    
    //action modes
    action_group_ = new QButtonGroup(tool_bar_);
    action_group_->setExclusive(true);
    
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(zoom));
    b->setToolTip("Action: Zoom + Translate");
    b->setShortcut(Qt::Key_Z);
    b->setCheckable(true);
    b->setWhatsThis("Action mode: Zoom + Translate<BR><BR>This mode allows to navigate in the data."
    								" The default is to zoom, Press the CTRL key for translation mode.<BR><BR>"
    								"A double-click resets the zoom.");
    action_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::AM_ZOOM);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(selection));
    b->setToolTip("Action: Select + Measure");
    b->setShortcut(Qt::Key_S);
    b->setCheckable(true);
    b->setWhatsThis("Action mode: Select + Measure<BR><BR>This mode allows to select peaks and"
    								" measure distances between peaks. The default is to select peaks. Press the"
    								" CTRL key for measurment mode.");
    action_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::AM_SELECT);
		tool_bar_->addWidget(b);

    connect(action_group_,SIGNAL(buttonClicked(int)),this,SLOT(setActionMode(int)));
    tool_bar_->addSeparator();
	   
    //intensity modes
    intensity_group_ = new QButtonGroup(tool_bar_);
    intensity_group_->setExclusive(true);
    
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(lin));
    b->setToolTip("Intensity: Normal");
    b->setShortcut(Qt::Key_N);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Normal<BR><BR>Intensity is displayed unmodified.");
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_NONE);
		tool_bar_->addWidget(b);
    
    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(logarithm));
    b->setToolTip("Intensity: Logarithmic");
    b->setShortcut(Qt::Key_L);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Logarithmic<BR><BR>Intensity is displayed in a logarithmic scale.");
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_LOG);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(percentage));
    b->setToolTip("Intensity: Percentage");
    b->setShortcut(Qt::Key_P);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Percentage<BR><BR>Intensity is displayed as a percentage of the layer"
    								" maximum intensity. If only one layer is displayed this mode behaves like the"
    								" normal mode. If more than one layer is displayed intensities are aligned.");
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_PERCENTAGE);
		tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QPixmap(snap));
    b->setToolTip("Intensity: Snap to maximum displayed intensity");
    b->setShortcut(Qt::Key_A);
    b->setCheckable(true);
    b->setWhatsThis("Intensity: Snap to maximum displayed intensity<BR><BR> In this mode the"
    								" color gradient is adapted to the maximum currently displayed intensity.");
    intensity_group_->addButton(b,SpectrumCanvas::SpectrumCanvas::IM_SNAP);
		tool_bar_->addWidget(b);

    connect(intensity_group_,SIGNAL(buttonClicked(int)),this,SLOT(setIntensityMode(int)));
    tool_bar_->addSeparator();

    //common buttons
    QAction* reset_zoom_button = tool_bar_->addAction(QPixmap(reset_zoom), "Reset Zoom", this, SLOT(resetZoom()));
    reset_zoom_button->setWhatsThis("Reset zoom: Zooms out as far as possible.");
    reset_zoom_button->setShortcut(Qt::Key_Backspace);

    tool_bar_->show();

		//--help toolbar--
		QToolBar* help_bar = addToolBar("Help tool bar");
		help_bar->addAction(QWhatsThis::createAction(help_bar));

    //--1D toolbar--
    tool_bar_1d_ = addToolBar("1D tool bar");

    //draw modes 1D
    draw_group_1d_ = new QButtonGroup(tool_bar_1d_);
    draw_group_1d_->setExclusive(true);
    
    b = new QToolButton(tool_bar_1d_);
    b->setIcon(QPixmap(peaks));
    b->setToolTip("Peak mode");
    b->setShortcut(Qt::Key_I);
    b->setCheckable(true);
    b->setWhatsThis("1D Draw mode: Peaks<BR><BR>Peaks are diplayed as sticks.");
    draw_group_1d_->addButton(b,Spectrum1DCanvas::DM_PEAKS);
		tool_bar_1d_->addWidget(b);
    
    b = new QToolButton(tool_bar_1d_);
    b->setIcon(QPixmap(lines));
    b->setToolTip("Raw data mode");
    b->setShortcut(Qt::Key_R);
    b->setCheckable(true);
    b->setWhatsThis("1D Draw mode: Raw data<BR><BR>Peaks are diplayed as a continous line.");
    draw_group_1d_->addButton(b,Spectrum1DCanvas::DM_CONNECTEDLINES);
		tool_bar_1d_->addWidget(b);

    connect(draw_group_1d_,SIGNAL(buttonClicked(int)),this,SLOT(setDrawMode1D(int)));
    tool_bar_->addSeparator();


    link_box_ = new QComboBox(tool_bar_1d_);
    link_box_->setToolTip("Linking spectra");
    link_box_->setWhatsThis("Linking spectra<BR><BR>Use this combobox to link two 1D spectra."
    												"Linked spectra zoom in/out together.");
    tool_bar_1d_->addWidget(link_box_);
    connect(link_box_,SIGNAL(activated(int)),this,SLOT(linkActiveTo(int)));

    //--2D toolbar--
    tool_bar_2d_ = addToolBar("2D tool bar");

    dm_surface_2d_ = tool_bar_2d_->addAction(QPixmap(colors),"Show colored surface");
    dm_surface_2d_->setCheckable(true);
    dm_surface_2d_->setWhatsThis("2D peak draw mode: Surface<BR><BR>The marching squares algorithm is applied"
    								             " to calculate a surface");
    connect(dm_surface_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    dm_contours_2d_ = tool_bar_2d_->addAction(QPixmap(contours),"Show contour lines");
    dm_contours_2d_->setCheckable(true);
    dm_contours_2d_->setWhatsThis("2D peak draw mode: Contour lines<BR><BR>The marching squares algorithm is applied"
    								              " to calculate contour lines.");
    connect(dm_contours_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    dm_precursors_2d_ = tool_bar_2d_->addAction(QPixmap(precursors),"Show MS/MS precursors");
    dm_precursors_2d_->setCheckable(true);
    dm_precursors_2d_->setWhatsThis("2D peak draw mode: Precursors<BR><BR>MS/MS precursor peaks are marked");
    connect(dm_precursors_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    projections_2d_ = tool_bar_2d_->addAction(QPixmap(projections), "Show Projections" ,this, SLOT(showProjections()));
    projections_2d_->setWhatsThis("Projections: Shows projections of peak data along RT and MZ axis.");


    dm_hull_2d_ = tool_bar_2d_->addAction(QPixmap(convexhulls),"Show feature convex hulls");
    dm_hull_2d_->setCheckable(true);
    dm_hull_2d_->setWhatsThis("2D feature draw mode: Convex hull<BR><BR>The convex hull of the feature is displayed");
    connect(dm_hull_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    dm_numbers_2d_ = tool_bar_2d_->addAction(QPixmap(numbers),"Show feature numbers");
    dm_numbers_2d_->setCheckable(true);
    dm_numbers_2d_->setWhatsThis("2D feature draw mode: Numbers<BR><BR>The feature number is displayed next to the feature");
    connect(dm_numbers_2d_, SIGNAL(toggled(bool)), this, SLOT(changeLayerFlag(bool)));

    //layer bar
    QDockWidget* layer_bar = new QDockWidget("Layers", this);
    addDockWidget(Qt::RightDockWidgetArea, layer_bar);
    layer_manager_ = new QListWidget(layer_bar);
    layer_manager_->setWhatsThis("Layer bar<BR><BR>Here the availabe layers are shown. You can select, hide"
    								              " and remove the layers using this bar.");

    layer_bar->setWidget(layer_manager_);
    layer_manager_->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(layer_manager_,SIGNAL(currentRowChanged(int)),this,SLOT(layerSelectionChange(int)));
		connect(layer_manager_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(layerContextMenu(const QPoint&)));
		connect(layer_manager_,SIGNAL(itemChanged(QListWidgetItem*)),this,SLOT(layerVisibilityChange(QListWidgetItem*)));
  }

  void TOPPViewBase::linkActiveTo(int index)
  {
  	Spectrum1DWidget* active = active1DWindow_();
  	if (active==0) return;
  	
  	//cout << "linkActiveTo() active: " << active->window_id << endl;
  	
  	//remove link if present
  	if (link_map_.find(active->window_id)!=link_map_.end())
  	{
  		SpectrumWidget* linked_to = window_(link_map_[active->window_id]);
  		if (linked_to != 0)
  		{
  			//cout << "  disconnect signals to: " << linked_to->window_id << endl;
	  		//disconnect outgoing signals
	  		disconnect(active->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),linked_to->canvas(),SLOT(setVisibleArea(DRange<2>)));
  			//disconnect incoming signals
  			disconnect(linked_to->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),active->canvas(),SLOT(setVisibleArea(DRange<2>)));
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
    SpectrumWidget* link_to = window_(link_box_->itemData(index).toInt());
		if (link_to!=0)
		{
			//cout << "  connecting signals to: " << link_to->window_id << endl;
		  connect(active->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),link_to->canvas(),SLOT(setVisibleArea(DRange<2>)));
      connect(link_to->canvas(),SIGNAL(visibleAreaChanged(DRange<2>)),active->canvas(),SLOT(setVisibleArea(DRange<2>)));
      //add entries to id map
      link_map_[link_to->window_id]=active->window_id;
      link_map_[active->window_id]=link_to->window_id;
		}
		else
		{
			cout << "linkActiveTo() - connect: Error, could not find window with id '" << link_box_->itemData(index).toInt() << "'!" << endl;
		}    
  }

  void TOPPViewBase::showStatusMessage(string msg, OpenMS::UInt time)
  {
    if (time==0)
    {
      message_label_->setText(msg.c_str());
      statusBar()->update();
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
    statusBar()->update();
  }

  void TOPPViewBase::changeGridLines()
  {
    SpectrumWidget* window = activeWindow_();
    if (window!=0)
    {
      window->canvas()->showGridLines(!window->canvas()->gridLinesShown());
    }
  }

  void TOPPViewBase::resetZoom()
  {
    SpectrumWidget* window = activeWindow_();
    if (window!=0)
    {
      window->canvas()->resetZoom();
    }
  }

  void TOPPViewBase::setActionMode(int index)
  {
    SpectrumWidget* w = activeWindow_();
    if (w)
    {
    	w->setActionMode((OpenMS::SpectrumCanvas::ActionModes)index);
  	}
  }

  void TOPPViewBase::setIntensityMode(int index)
  {
    SpectrumWidget* w = activeWindow_();
    if (w)
    {
    	w->setIntensityMode((OpenMS::SpectrumCanvas::IntensityModes)index);
  	}
  }

  void TOPPViewBase::setDrawMode1D(int index)
  {
    Spectrum1DWidget* w = active1DWindow_();
    if (w)
    {
    	w->canvas()->setDrawMode((OpenMS::Spectrum1DCanvas::DrawModes)index);
  	}
  }

  void TOPPViewBase::changeLayerFlag(bool on)
  {
		QAction* action = qobject_cast<QAction *>(sender());
    if (Spectrum2DWidget* win = active2DWindow_())
    {
    	//peaks
	    if (action == dm_contours_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::P_CONTOURS,on);
			}
			else if (action == dm_surface_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::P_SURFACE,on);
			}
			else if (action == dm_precursors_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::P_PRECURSORS,on);
			}
			//features
			else if (action == dm_hull_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::F_HULLS,on);
			}
			else if (action == dm_numbers_2d_)
			{
		    win->canvas()->setLayerFlag(LayerData::F_NUMBERS,on);
			}
		}
  }

  void TOPPViewBase::updateToolbar()
  {
    SpectrumWidget* w = activeWindow_();

    if (w)
    {
      //set action mode
     	action_group_->button(w->getActionMode())->setChecked(true);

      //set intensity mode
     	intensity_group_->button(w->canvas()->getIntensityMode())->setChecked(true);
    }

    //1D
    Spectrum1DWidget* w1 = active1DWindow_();
    if (w1)
    {
      //draw mode
      draw_group_1d_->button(w1->canvas()->getDrawMode())->setChecked(true);

      //update link selector
      int item_index = -1;
      link_box_->clear();
      link_box_->insertItem(++item_index,"<unlinked>",0);
      QWidgetList windows = ws_->windowList();
      for ( int i = 0; i < windows.count(); ++i )
      {
        Spectrum1DWidget* window = dynamic_cast<Spectrum1DWidget*>(windows.at(i));
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
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_SELECT)->setEnabled(true);
    }

    //2d
    Spectrum2DWidget* w2 = active2DWindow_();
    if (w2)
    {
      //peak draw modes
      if (w2->canvas()->getCurrentLayer().type == LayerData::DT_PEAK)
      {
      	dm_surface_2d_->setVisible(true);
      	dm_contours_2d_->setVisible(true);
      	dm_precursors_2d_->setVisible(true);
      	projections_2d_->setVisible(true);
      	dm_hull_2d_->setVisible(false);
      	dm_numbers_2d_->setVisible(false);
      	dm_surface_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::P_SURFACE));
      	dm_contours_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::P_CONTOURS));
				dm_precursors_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::P_PRECURSORS));
			}
			//feature draw modes
			else
			{
      	dm_surface_2d_->setVisible(false);
      	dm_contours_2d_->setVisible(false);
      	dm_precursors_2d_->setVisible(false);
      	projections_2d_->setVisible(false);
      	dm_hull_2d_->setVisible(true);
      	dm_numbers_2d_->setVisible(true);
      	dm_hull_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::F_HULLS));
      	dm_numbers_2d_->setChecked(w2->canvas()->getLayerFlag(LayerData::F_NUMBERS));
			}
      //show/hide toolbars and buttons
      tool_bar_1d_->hide();
      tool_bar_2d_->show();
      action_group_->button(SpectrumCanvas::SpectrumCanvas::AM_SELECT)->setEnabled(true);
    }

    //1D
    Spectrum3DWidget* w3 = active3DWindow_();
    if (w3)
    {
      //show/hide toolbars and buttons
      tool_bar_1d_->hide();
      tool_bar_2d_->hide();
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
    for (UInt i = 0; i<cc->getLayerCount(); ++i)
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
			context_menu->addAction("Delete");
			if (context_menu->exec(layer_manager_->mapToGlobal(pos)))
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
  		tab_bar_->setCurrentIndex(dynamic_cast<SpectrumWidget*>(w)->window_id);
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

  void TOPPViewBase::showAsWindow_(SpectrumWidget* sw, const String& caption)
  {
  	ws_->addWindow(sw);
    connect(sw->canvas(),SIGNAL(layerActivated(QWidget*)),this,SLOT(updateToolbar()));
    connect(sw,SIGNAL(sendStatusMessage(std::string,OpenMS::UInt)),this,SLOT(showStatusMessage(std::string,OpenMS::UInt)));
    connect(sw,SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
    connect(sw,SIGNAL(modesChanged(QWidget*)),this,SLOT(updateToolbar()));
  
  	Spectrum2DWidget* sw2 = dynamic_cast<Spectrum2DWidget*>(sw);
  	if (sw2 != 0)
  	{
  		connect(sw2->getHorizontalProjection(),SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
  		connect(sw2->getVerticalProjection(),SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
  		connect(sw2,SIGNAL(showCurrentPeaksAs3D()),this,SLOT(showCurrentPeaksAs3D()));
  		connect(sw2,SIGNAL(showSpectrumAs1D(int)),this,SLOT(showSpectrumAs1D(int)));
  	}
  	
	  sw->setWindowTitle(caption.c_str());
    addTab_(sw,caption);
  }

  void TOPPViewBase::gotoDialog()
  {
    SpectrumWidget* w = activeWindow_();
    if (w)
    {
      w->showGoToDialog();
    }
  }

  SpectrumWidget*  TOPPViewBase::activeWindow_() const
  {
    return dynamic_cast<SpectrumWidget*>(ws_->activeWindow());
  }

  SpectrumCanvas*  TOPPViewBase::activeCanvas_() const
  {
    SpectrumWidget* sw = dynamic_cast<SpectrumWidget*>(ws_->activeWindow());
    if (sw == 0)
    {
    	return 0;
    }
    return sw->canvas();
  }

  Spectrum1DWidget* TOPPViewBase::active1DWindow_() const
  {
    Spectrum1DWidget* s1;
    if ((s1 = dynamic_cast<Spectrum1DWidget*>(ws_->activeWindow())))
    {
      return s1;
    }
    return 0;
  }

  Spectrum2DWidget* TOPPViewBase::active2DWindow_() const
  {
    Spectrum2DWidget* s2;
    if ((s2 = dynamic_cast<Spectrum2DWidget*>(ws_->activeWindow())))
    {
      return s2;
    }
    return 0;
  }

  Spectrum3DWidget* TOPPViewBase::active3DWindow_() const
  {
    Spectrum3DWidget* s3;
    if ((s3 = dynamic_cast<Spectrum3DWidget*>(ws_->activeWindow())))
    {
      return s3;
    }
    return 0;
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

    //load preferences, if file exists
    if (File::exists(filename))
    {
    	Param tmp;
    	tmp.load(filename);
    	//apply preferences if they are of the current TOPPView version
    	if(tmp.getValue("Preferences:Version").toString()==VersionInfo::getVersion())
    	{
      	setParameters(tmp);
    	}
    	else
    	{
    		cout << "The preferences files '" << filename  
    		     << "' is replaced as it is not of the current TOPPView version." << endl;
    	}
    }
    else
    {
      if (filename != default_ini_file)
      {
        cerr << "Unable to load INI File: '" << filename << "'" << endl;
      }
    }
    param_.setValue("PreferencesFile" , filename);

    //set the recent files
    Param p = param_.copy("Preferences:RecentFiles");
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
    param_.remove("Preferences:RecentFiles");

    for (int i = 0; i < recent_files_.size(); ++i)
    {
      param_.setValue("Preferences:RecentFiles:"+String(i),recent_files_[i].toStdString());
    }
		
		//set version
		param_.setValue("Preferences:Version",VersionInfo::getVersion());
		
    //save only the subsection that begins with "Preferences:"
    try
    {
      param_.copy("Preferences:").store(string(param_.getValue("PreferencesFile")));
    }
    catch(Exception::UnableToCreateFile& e)
    {
      cerr << "Unable to create INI File: '" << string(param_.getValue("PreferencesFile")) << "'" << endl;
    }
  }

  void TOPPViewBase::openRecentFile()
  {
		QAction* action = qobject_cast<QAction *>(sender());
    if (action)
		{
	  	setCursor(Qt::WaitCursor);
 	  	OpenDialog::Mower mow = OpenDialog::NO_MOWER;
			if ( (String)param_.getValue("Preferences:MapIntensityCutoff")=="Noise Estimator")
			{
				mow = OpenDialog::NOISE_ESTIMATOR;
			}
   		addSpectrum(action->text().toStdString(),!recent_as_new_layer_->isChecked(),(String)param_.getValue("Preferences:DefaultMapView")=="2D",true,mow);
			setCursor(Qt::ArrowCursor); 	
		}
 }

  void TOPPViewBase::openSpectrumDialog()
  {
    OpenDialog dialog(param_,this);
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

	void TOPPViewBase::editParamDialog()
	{
		// CREATE DIALOG
		QDialog dialog(this);
		QGridLayout* layout = new QGridLayout(&dialog);
		//Editor
		ParamEditor* edit = new ParamEditor(&dialog);		
		layout->addWidget(edit,0,0,1,3);
		//Stretch
		layout->setColumnStretch(0,2);
		//Store button
		QPushButton* button = new QPushButton("Cancel",&dialog);
		connect(button,SIGNAL(pressed()),&dialog,SLOT(reject()));
		layout->addWidget(button,1,1);
		//Cancel button
		button = new QPushButton("OK",&dialog);
		connect(button,SIGNAL(pressed()),&dialog,SLOT(accept()));
		layout->addWidget(button,1,2);
		
		//LOAD DATA	
		QString name = QFileDialog::getOpenFileName(this,"Select a INI file",param_.getValue("Preferences:DefaultPath").toString().c_str(),tr("ini files (*.ini);; all files (*.*)"));
		if (name=="") return;

		Param p;
		p.load(name.toAscii().data());
		edit->loadEditable(p);

		//EXECUTE DIALOG + STORE DATA
		if (dialog.exec())
		{
			if (!File::writable(name.toAscii().data()))
			{
				QMessageBox::critical(this,"Error writing file!",(String("Cannot write to '")+name.toAscii().data()+"'!").c_str());
				return;
			}
			edit->store();
			p.store(name.toAscii().data());
		}
	}

	void TOPPViewBase::showTOPPDialog()
	{
		String tmp_dir = param_.getValue("Preferences:TmpPath").toString();
		
		ToolsDialog dialog(this,tmp_dir);
	
		if(dialog.exec()==QDialog::Accepted)
		{
			//check if there is a active window
			if (ws_->activeWindow())
			{
				const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
				//warn if hidden layer => wrong layer selected...
				if (!layer.visible)
				{
					QMessageBox::warning(this,"Warning","The current layer is not visible!");
				}
				if (!File::writable(tmp_dir+"/in"))
				{
					QMessageBox::critical(this,"Error creating temporary file!",(String("Cannot write to '")+tmp_dir+"/in'!").c_str());
					return;
				}
				if (!File::writable(tmp_dir+"/out"))
				{
					QMessageBox::critical(this,"Error creating temporary file!",(String("Cannot write to '")+tmp_dir+"/out'!").c_str());
					return;
				}
				if (layer.type==LayerData::DT_PEAK)
				{
					MzDataFile().store(tmp_dir+"/in",layer.peaks);
				}
				else if (layer.type==LayerData::DT_FEATURE)
				{
					FeatureMapFile().store(tmp_dir+"/in",layer.features);
				}
				else if (layer.type==LayerData::DT_FEATURE_PAIR)
				{
					LayerData::FeatureMapType feature_map=layer.features;
					FeatureMap<>::ConstIterator end=--feature_map.end();
					vector< ElementPair< Feature > > feature_pairs;
					FeatureMap<>::ConstIterator next;
					for(FeatureMap<>::ConstIterator i=feature_map.begin();i<end;i+=2)
					{
						next=++i;
						--i;
						feature_pairs.push_back(ElementPair<>(*i,*next));
					}
					FeaturePairsFile().store(tmp_dir+"/in",feature_pairs);
				}
				else
				{
					return;
				}
					
				String output = dialog.isOutputOnly() ? String(""):String("-"+dialog.getOutput()+" "+tmp_dir+"/out");
				String call = dialog.getTool() + " -ini "+tmp_dir+"/in.ini -"+dialog.getInput()+" "+tmp_dir+"/in "+output +" > "+tmp_dir+"/TOPP.log 2>&1";
					
				if (system(call.c_str())!=0)
				{
					TextFile f;
					f.load(tmp_dir+"/TOPP.log");
					String log_file;
					log_file.implode(f.begin(),f.end(),"<BR>");
					QMessageBox::critical(this,"Execution of TOPP tool not successful!",log_file.c_str());
				}
				else if (!File::readable(tmp_dir+"/out"))
				{
					QMessageBox::critical(this,"Error creating temporary file!",(String("Cannot read '")+tmp_dir+"/in'!").c_str());
					return;
				}
				else if(dialog.isWindow())
				{
					addSpectrum(tmp_dir+"/out",true,true,true);
				}
				else if(dialog.isLayer())
				{
					addSpectrum(tmp_dir+"/out",false,true,true);
				}
				else if(dialog.isOutputOnly())
				{
					TextFile f;
					f.load(tmp_dir+"/TOPP.log");
					String log_file;
					log_file.implode(f.begin(),f.end(),"<BR>");
					QMessageBox::information(this,("Standard Output of Tool: "+dialog.getTool()).c_str(), log_file.c_str());
				}
			}
		}
	}
	
	const LayerData* TOPPViewBase::getCurrentLayer() const
	{
		SpectrumCanvas* canvas = activeCanvas_();
		if (canvas==0)
		{
			return 0;
		}
		return &(canvas->getCurrentLayer());
	}

  void TOPPViewBase::showProjections()
  {
    Spectrum2DWidget* w = active2DWindow_();
    if (w)
    {
    	w->canvas()->showProjections();
    }
  }

	void TOPPViewBase::annotateWithID()
	{
		QString name = QFileDialog::getOpenFileName(this,"Select identification data",param_.getValue("Preferences:DefaultPath").toString().c_str(),tr("identfication files (*.analysisXML);; all files (*.*)"));
		
		if(name!="")
		{
			//check if there is a active window
			if (ws_->activeWindow())
			{
				const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
				//warn if hidden layer => wrong layer selected...
				if (!layer.visible)
				{
					QMessageBox::warning(this,"Warning","The current layer is not visible!");
				}
				
				//load id data
				vector<IdentificationData> identifications; 
				vector<ProteinIdentification> protein_identifications; 
				AnalysisXMLFile().load(name.toStdString(), protein_identifications, identifications);
				if (layer.type==LayerData::DT_PEAK)
				{
					IDSpectrumMapper().annotate(const_cast<LayerData&>(layer).peaks, identifications, 0.1);
				}
				else if (layer.type==LayerData::DT_FEATURE || layer.type==LayerData::DT_FEATURE_PAIR)
				{
					IDFeatureMapper().annotate(const_cast<LayerData&>(layer).features,identifications,protein_identifications);
				}
			}
		}
	}

	void TOPPViewBase::showCurrentPeaksAs3D()
	{
  	//check if there is a active window
    if (ws_->activeWindow())
    {
      const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
    	const SpectrumCanvas::AreaType& area = activeWindow_()->canvas()->getVisibleArea();
    	
    	if (layer.type==LayerData::DT_PEAK)
    	{
    		//open new 3D widget
    		Spectrum3DWidget* w = new Spectrum3DWidget(param_, ws_);
  			SpectrumCanvas::ExperimentType& out = w->canvas()->addEmptyPeakLayer();
  			
    		for (LayerData::ExperimentType::ConstIterator it=layer.peaks.RTBegin(area.min()[1]); it!=layer.peaks.RTBegin(area.max()[1]); ++it)
    		{
    			if (it->getMSLevel()!=1) continue;
    			SpectrumCanvas::ExperimentType::SpectrumType spectrum;
    				
  				spectrum.setRT(it->getRT());
  				spectrum.setMSLevel(it->getMSLevel());
  				for (LayerData::ExperimentType::SpectrumType::ConstIterator it2 = it->MZBegin(area.min()[0]); it2!= it->MZEnd(area.max()[0]); ++it2)
  				{
  					if ( it2->getIntensity() >= layer.min_int && it2->getIntensity() <= layer.max_int)
  					{
  						spectrum.push_back(*it2);
  					}
  				}
  				out.push_back(spectrum);
    		}
    		
    		if (out.size()==0) //no extracted data
    		{
    		  QMessageBox::critical(this,"Error","The displayed region of the current layer is empty!");
    		  delete(w);
    		}
    		else //finish adding
    		{
    			String caption = layer.name + " (3D)";
    			w->canvas()->finishAdding(0.0);
					w->canvas()->setCurrentLayerName(caption);
		      showAsWindow_(w,caption);
		      w->showMaximized();
			    updateLayerbar();
    		}
    		
			}
    }
	}

	void TOPPViewBase::showSpectrumAs1D(int index)
	{
  	//check if there is a active window
    if (ws_->activeWindow())
    {
      const LayerData& layer = activeWindow_()->canvas()->getCurrentLayer();
    	
    	if (layer.type==LayerData::DT_PEAK)
    	{
    		//open new 1D widget
    		Spectrum1DWidget* w = new Spectrum1DWidget(param_, ws_);
  			w->canvas()->addEmptyPeakLayer().push_back(layer.peaks[index]);
  			String caption = layer.name + " (RT: " + layer.peaks[index].getRT() + ")";
  			w->canvas()->finishAdding(0.0);
				w->canvas()->setCurrentLayerName(caption);
	      showAsWindow_(w,caption);
	      w->showMaximized();
		    updateLayerbar();
			}
    }
	}

} //namespace OpenMS


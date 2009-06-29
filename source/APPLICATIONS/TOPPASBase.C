// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Johannes Junker $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPASBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h> 
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASWidget.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>

//Qt
#include <QtGui/QToolBar>
#include <QtGui/QDockWidget>
#include <QtGui/QListWidget>
#include <QtGui/QListWidgetItem>
#include <QtGui/QTreeWidget>
#include <QtGui/QTreeWidgetItem>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QToolButton>
#include <QtGui/QMessageBox>
#include <QtGui/QToolTip>
#include <QtGui/QFileDialog>
#include <QtGui/QWhatsThis>
#include <QtGui/QInputDialog>
#include <QtGui/QTextEdit>
#include <QtGui/QCheckBox>
#include <QtGui/QCloseEvent>
#include <QtGui/QDesktopServices>
#include <QtCore/QUrl>
#include <QtGui/QSplashScreen>
#include <QtGui/QVBoxLayout>
#include <QtGui/QApplication>
#include <QtGui/QLabel>

using namespace std;

namespace OpenMS
{
  using namespace Internal;

  TOPPASBase::TOPPASBase(QWidget* parent):
      QMainWindow(parent),
      DefaultParamHandler("TOPPASBase")
  {
  	setWindowTitle("TOPPAS");
    //setWindowIcon(QIcon(toppas)); TODO
    //prevents errors caused by too small width,height values
    setMinimumSize(400,400);

    // create dummy widget (to be able to have a layout), Tab bar and workspace
    QWidget* dummy = new QWidget(this);
    setCentralWidget(dummy);
    QVBoxLayout* box_layout = new QVBoxLayout(dummy);
    tab_bar_ = new EnhancedTabBar(dummy);
    tab_bar_->setWhatsThis("Tab bar<BR><BR>Close tabs through the context menu or by double-clicking them.<BR>The tab bar accepts drag-and-drop from the layer bar.");
    tab_bar_->addTab("dummy",1336);
    tab_bar_->setMinimumSize(tab_bar_->sizeHint());
    tab_bar_->removeId(1336);
    //connect slots and sigals for selecting spectra
    connect(tab_bar_,SIGNAL(currentIdChanged(int)),this,SLOT(focusByTab(int)));
    connect(tab_bar_,SIGNAL(aboutToCloseId(int)),this,SLOT(closeByTab(int)));

    box_layout->addWidget(tab_bar_);
    ws_=new QWorkspace(dummy);
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateTabBar(QWidget*)));
    connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateMenu()));

    box_layout->addWidget(ws_);
		
		//################## MENUS #################
    // File menu
    QMenu* file = new QMenu("&File",this);
    menuBar()->addMenu(file);
    file->addAction("&New",this,SLOT(newFileDialog()), Qt::CTRL+Qt::Key_N);
    file->addAction("&Open",this,SLOT(openFileDialog()), Qt::CTRL+Qt::Key_O);
    file->addAction("&Save",this,SLOT(saveFileDialog()), Qt::CTRL+Qt::Key_S);
		file->addAction("Save &As",this,SLOT(saveFileDialog()), Qt::CTRL+Qt::SHIFT+Qt::Key_S);
    file->addAction("&Close",this,SLOT(closeFile()), Qt::CTRL+Qt::Key_W);
		file->addSeparator();

    file->addSeparator();
    file->addAction("&Preferences",this, SLOT(preferencesDialog()));
    file->addAction("&Quit",qApp,SLOT(quit()));

    //Advanced menu
    //QMenu* advanced = new QMenu("&Advanced",this);
    //menuBar()->addMenu(advanced);
    //advanced->addAction("&Refresh definitions",this,SLOT(refreshDefinitions()), Qt::CTRL+Qt::Key_R);

		//Windows menu
    QMenu * windows = new QMenu("&Windows", this);
    menuBar()->addMenu(windows);

		//Help menu
		QMenu* help = new QMenu("&Help", this);
		menuBar()->addMenu(help);
		//help->addAction(QWhatsThis::createAction(help));
		//help->addSeparator();
		QAction* action = help->addAction("OpenMS website",this,SLOT(showURL()));
		action->setData("http://www.OpenMS.de");
		help->addSeparator();
		help->addAction("&About",this,SLOT(showAboutDialog()));

    //create status bar
    message_label_ = new QLabel(statusBar());
    statusBar()->addWidget(message_label_,1);

		//################## TOOLBARS #################
    //create toolbars and connect signals

  	//--Basic tool bar--
    //tool_bar_ = addToolBar("Basic tool bar");
    //tool_bar_->show();

		//################## Dock widgets #################
    
    //TOPP tools window
    QDockWidget* topp_tools_bar = new QDockWidget("TOPP", this);
    addDockWidget(Qt::LeftDockWidgetArea, topp_tools_bar);
    tools_tree_view_ = new QTreeWidget(topp_tools_bar);
		connect (tools_tree_view_, SIGNAL(itemPressed(QTreeWidgetItem*, int)), this, SLOT(treeViewItemPressed_(QTreeWidgetItem*, int)));
    tools_tree_view_->setWhatsThis("TOPP tools list<BR><BR>All available TOPP tools are shown here.");
    tools_tree_view_->setColumnCount(1);
  	QStringList header_labels;
  	header_labels.append(QString("Tool"));
  	tools_tree_view_->setHeaderLabels(header_labels);
    topp_tools_bar->setWidget(tools_tree_view_);
    
    QTreeWidgetItem* item = new QTreeWidgetItem((QTreeWidget*)0);
    item->setText(0, "<Input file>");
    tools_tree_view_->addTopLevelItem(item);
    item = new QTreeWidgetItem((QTreeWidget*)0);
    item->setText(0, "<Input file list>");
    tools_tree_view_->addTopLevelItem(item);
    item = new QTreeWidgetItem((QTreeWidget*)0);
    item->setText(0, "<Output file>");
    tools_tree_view_->addTopLevelItem(item);
    item = new QTreeWidgetItem((QTreeWidget*)0);
    item->setText(0, "<Output file list>");
    tools_tree_view_->addTopLevelItem(item);
    
    
    Map<String,StringList> tools_list = TOPPBase::getToolList();
    QTreeWidgetItem* parent_item;
    for (Map<String,StringList>::Iterator it = tools_list.begin(); it != tools_list.end(); ++it)
    {
    	item = new QTreeWidgetItem((QTreeWidget*)0);
    	item->setText(0, it->first.toQString());
    	tools_tree_view_->addTopLevelItem(item);
   		parent_item = item;
   		StringList types = it->second;
   		for (StringList::iterator types_it = types.begin(); types_it != types.end(); ++types_it)
   		{
   			item = new QTreeWidgetItem(parent_item);
   			item->setText(0, types_it->toQString());
   		}
    }
    tools_tree_view_->resizeColumnToContents(0);
    
		windows->addAction(topp_tools_bar->toggleViewAction());
		
    //Blocks window
//     QDockWidget* blocks_bar = new QDockWidget("Blocks", this);
//     addDockWidget(Qt::LeftDockWidgetArea, blocks_bar);
//     blocks_list_ = new QListWidget(blocks_bar);
//     blocks_list_->setWhatsThis("Blocks list<BR><BR>Custom analysis pipelines are shown here. They can be used as if they were TOPP tools themselves.");
//     blocks_bar->setWidget(blocks_list_);

		//log window
		QDockWidget* log_bar = new QDockWidget("Log", this);
		addDockWidget(Qt::BottomDockWidgetArea, log_bar);
		log_ = new QTextEdit(log_bar);
		log_->setReadOnly(true);
		log_bar->setWidget(log_);
		log_bar->hide();
    //windows->addAction("&Show log window",log_bar,SLOT(show()));
		windows->addAction(log_bar->toggleViewAction());
		//################## DEFAULTS #################
    //general
   	defaults_.setValue("preferences:default_path", ".", "Default path for loading and storing files.");
    defaults_.setValue("preferences:default_path_current", "true", "If the current path is preferred over the default path.");
		defaults_.setValidStrings("preferences:default_path_current",StringList::create("true,false"));
		defaults_.setValue("preferences:tmp_file_path", QDir::tempPath(), "Path where temporary files can be created.");
		defaults_.setValue("preferences:version","none","OpenMS version, used to check if the TOPPAS.ini is up-to-date");
		
  	defaultsToParam_();

  	//load param file
    loadPreferences();
		
		//set current path
		current_path_ = param_.getValue("preferences:default_path");
		
  	//update the menu
  	updateMenu();
	}

  TOPPASBase::~TOPPASBase()
  {
  	savePreferences();
  }

	void TOPPASBase::refreshDefinitions()
	{
	
	}
	
	void TOPPASBase::openFileDialog()
  {
		
  }
  
  void TOPPASBase::newFileDialog()
  {
  	TOPPASWidget* tw = new TOPPASWidget(Param(), ws_);
  	showAsWindow_(tw, "New");
  }
	
	void TOPPASBase::saveFileDialog()
	{
		
	}
	
	void TOPPASBase::preferencesDialog()
  {
		// do something...
		savePreferences();
  }
	
	void TOPPASBase::showAsWindow_(TOPPASWidget* tw, const String& caption)
  {
  	ws_->addWindow(tw);
    connect(tw,SIGNAL(sendStatusMessage(std::string,OpenMS::UInt)),this,SLOT(showStatusMessage(std::string,OpenMS::UInt)));
    connect(tw,SIGNAL(sendCursorStatus(double,double)),this,SLOT(showCursorStatus(double,double)));
		connect(tw,SIGNAL(toolDroppedOnWidget(double,double)),this,SLOT(insertNewVertex_(double,double)));
	  tw->setWindowTitle(caption.toQString());

		//add tab with id
  	static int window_counter = 1337;
  	tw->window_id = window_counter++;

    tab_bar_->addTab(caption.toQString(), tw->window_id);

    //connect slots and sigals for removing the widget from the bar, when it is closed
    //- through the menu entry
    //- through the tab bar
    //- thourgh the MDI close button
    connect(tw,SIGNAL(aboutToBeDestroyed(int)),tab_bar_,SLOT(removeId(int)));

    tab_bar_->setCurrentId(tw->window_id);

		//show first window maximized (only visible windows are in the list)
		if (ws_->windowList().count()==0)
		{
			tw->showMaximized();
		}
		else
		{
			tw->show();
		}
  }

	
  void TOPPASBase::closeEvent(QCloseEvent* event)
  {
  	ws_->closeAllWindows();
  	event->accept();
  }

	void TOPPASBase::showURL()
	{
		QAction* action = qobject_cast<QAction*>(sender());
		if (!QDesktopServices::openUrl(QUrl(action->data().toString())))
		{
			QMessageBox::warning(this, tr("Error"),
												 tr("Unable to open\n") +
														action->data().toString() +
														tr("\n\nPossible reason: security settings or misconfigured Operating System"));
		}
	}
	
  TOPPASWidget* TOPPASBase::window_(int id) const
  {
		//cout << "Looking for tab with id: " << id << endl;
  	QList<QWidget*> windows = ws_->windowList();
		for(int i=0; i< windows.size(); ++i)
		{
			TOPPASWidget* window = qobject_cast<TOPPASWidget*>(windows.at(i));
			//cout << "  Tab " << i << ": " << window->window_id << endl;
			if (window->window_id == id)
			{
				return window;
			}
		}
		return 0;
  }
  
  TOPPASWidget* TOPPASBase::activeWindow_() const
  {
  	if (!ws_->activeWindow()) return 0;
    return qobject_cast<TOPPASWidget*>(ws_->activeWindow());
  }

  void TOPPASBase::closeByTab(int id)
  {
  	TOPPASWidget* window = window_(id);
  	if (window)
  	{
  		window->close();
  		updateMenu();
  	}
  }

  void TOPPASBase::focusByTab(int id)
  {
  	TOPPASWidget* window = window_(id);
  	if (window)
  	{
  		window->setFocus();
  	}
  }

  void TOPPASBase::closeFile()
  {
    ws_->activeWindow()->close();
    updateMenu();
  }

  void TOPPASBase::showStatusMessage(string msg, OpenMS::UInt time)
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

	void TOPPASBase::showCursorStatus(double /*x*/, double /*y*/)
	{
		// TODO
	}

  void TOPPASBase::updateToolBar()
  {
    
  }

  void TOPPASBase::updateTabBar(QWidget* w)
  {
  	if (w)
  	{
  		Int window_id = qobject_cast<TOPPASWidget*>(w)->window_id;
  		tab_bar_->setCurrentId(window_id);
  	}
  }

  void TOPPASBase::loadPreferences(String filename)
  {
    //compose default ini file path
    String default_ini_file = String(QDir::homePath()) + "/.TOPPAS.ini";

    if (filename=="")
    {
      filename = default_ini_file;
    }

    //load preferences, if file exists
    if (File::exists(filename))
    {
    	bool error = false;
    	Param tmp;
    	tmp.load(filename);
    	//apply preferences if they are of the current TOPPAS version
    	if(tmp.exists("preferences:version") && tmp.getValue("preferences:version").toString()==VersionInfo::getVersion())
    	{
    		try
    		{
      		setParameters(tmp);
    		}
    		catch (Exception::InvalidParameter& /*e*/)
    		{
    			error = true;
    		}
    	}
    	else
    	{
				error = true;
    	}
			//set parameters to defaults when something is fishy with the parameters file
			if (error)
			{
  			//reset parameters
  			setParameters(Param());

				cerr << "The TOPPAS preferences files '" << filename << "' was ignored. It is no longer compatible with this TOPPAS version and will be replaced." << endl;
			}
    }
    else if (filename != default_ini_file)
    {
    	cerr << "Unable to load INI File: '" << filename << "'" << endl;
    }
    param_.setValue("PreferencesFile" , filename);
  }

  void TOPPASBase::savePreferences()
  {
		//set version
		param_.setValue("preferences:version",VersionInfo::getVersion());

    //save only the subsection that begins with "preferences:"
    try
    {
      param_.copy("preferences:").store(string(param_.getValue("PreferencesFile")));
    }
    catch(Exception::UnableToCreateFile& /*e*/)
    {
      cerr << "Unable to create INI File: '" << string(param_.getValue("PreferencesFile")) << "'" << endl;
    }
  }

	void TOPPASBase::showAboutDialog()
	{
		//dialog and grid layout
		QDialog* dlg = new QDialog(this);
		QGridLayout* grid = new QGridLayout(dlg);
		dlg->setWindowTitle("About TOPPAS");

		//image TODO
		//QLabel* label = new QLabel(dlg);
		//QPixmap image(Oesterberg);
		//label->setPixmap(image);
		//grid->addWidget(label,0,0);

		//text
		QString text = QString("<BR>"
									 				 "<FONT size=+3>TOPPAS</font><BR>"
									 				 "<BR>"
													 "Version: %1<BR>"
													 "<BR>"
													 "OpenMS and TOPP is free software available under the<BR>"
													 "Lesser GNU Public License (LGPL)<BR>"
													 "<BR>"
													 "<BR>"
													 "<BR>"
													 "<BR>"
													 "<BR>"
													 "Any published work based on TOPP and OpenMS shall cite these papers:<BR>"
													 "Sturm et al., BMC Bioinformatics (2008), 9, 163<BR>"
													 "Kohlbacher et al., Bioinformatics (2007), 23:e191-e197<BR>"
													 ).arg(VersionInfo::getVersion().toQString());
		QLabel* label = new QLabel(text,dlg);
		grid->addWidget(label,0,1,Qt::AlignTop | Qt::AlignLeft);

		//execute
		dlg->exec();
	}
	
  void TOPPASBase::updateMenu()
  {
  	
  }

  void TOPPASBase::showLogMessage_(TOPPASBase::LogState state, const String& heading, const String& body)
  {
		//Compose current time string
		DateTime d;
		d.now();

		String state_string;
		switch(state)
		{
			case LS_NOTICE: state_string = "NOTICE"; break;
			case LS_WARNING: state_string = "WARNING"; break;
			case LS_ERROR: state_string = "ERROR"; break;
		}

		//update log
		log_->append("==============================================================================");
		log_->append((d.getTime() + " " + state_string + ": " + heading).toQString());
		log_->append(body.toQString());

  	//show log tool window
		qobject_cast<QWidget *>(log_->parent())->show();
  }

	void TOPPASBase::keyPressEvent(QKeyEvent* /*e*/)
	{
 		
	}
	
	void TOPPASBase::updateCurrentPath()
	{
		//do not update if the user disabled this feature.
		if (param_.getValue("preferences:default_path_current")!="true") return;
		
		//reset
		current_path_ = param_.getValue("preferences:default_path");
		
		//update if the current layer has a path associated TODO
		//if (activeCanvas_() && activeCanvas_()->getLayerCount()!=0 && activeCanvas_()->getCurrentLayer().filename!="")
		//{
		//	current_path_ = File::path(activeCanvas_()->getCurrentLayer().filename);
		//}
	}
	
	void TOPPASBase::insertNewVertex_(double x, double y)
	{
		QTreeWidgetItem* current_tool = tools_tree_view_->currentItem();
		String tool_name = String(current_tool->text(0));
		TOPPASVertex* tv = 0;
		
		if (tool_name == "<Input file list>")
		{
			tv = new TOPPASInputFileListVertex();
		}
		else if (tool_name == "<Input file>")
		{
			tv = new TOPPASInputFileVertex();
		}
		else if (tool_name == "<Output file list>")
		{
			tv = new TOPPASOutputFileListVertex();
		}
		else if (tool_name == "<Output file>")
		{
			tv = new TOPPASOutputFileVertex();
		}
		else // node is a TOPP tool
		{	
			if (current_tool->childCount() > 0)
			{
				// tool has a type, but is selected itself (instead of type)
				return;
			}
			String tool_type;
			if (current_tool->parent() != 0)
			{
				// selected item is a type
				tool_type = String(current_tool->text(0));
				tool_name = String(current_tool->parent()->text(0));
			}
			else
			{
				// normal tool which does not have type selected
				tool_name = String(current_tool->text(0));
				tool_type = "";
			}
			
			tv = new TOPPASToolVertex(tool_name, tool_type);
		}
		
		tv->setPos(x,y);
		activeWindow_()->getScene()->addVertex(tv);
		
		connect(tv,SIGNAL(clicked()),activeWindow_()->getScene(),SLOT(itemClicked()));
		connect(tv,SIGNAL(doubleClicked()),activeWindow_()->getScene(),SLOT(itemDoubleClicked()));
		connect(tv,SIGNAL(hoveringEdgePosChanged(const QPointF&)),activeWindow_()->getScene(),SLOT(updateHoveringEdgePos(const QPointF&)));
		connect(tv,SIGNAL(newHoveringEdge(const QPointF&)),activeWindow_()->getScene(),SLOT(addHoveringEdge(const QPointF&)));
		connect(tv,SIGNAL(finishHoveringEdge()),activeWindow_()->getScene(),SLOT(finishHoveringEdge()));
	}
	
	void TOPPASBase::treeViewItemPressed_(QTreeWidgetItem* item, int /*column*/)
	{
		if (item == 0)
		{
			return;
		}
		
		if (item->childCount() > 0)
		{
			// a tool with types was pressed but only its types are allowed to be dragged
			tools_tree_view_->setDragEnabled(false);
		}
		else
		{
			tools_tree_view_->setDragEnabled(true);
		}
	}


} //namespace OpenMS


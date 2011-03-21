// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/VISUAL/APPLICATIONS/TOPPASBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h> 
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASWidget.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASTabBar.h>
#include <OpenMS/VISUAL/TOPPASResources.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

//Qt
#include <QtGui/QToolBar>
#include <QtGui/QDesktopWidget>
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
#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QSet>
#include <QtCore/QMap>

using namespace std;

namespace OpenMS
{
  using namespace Internal;

	int TOPPASBase::node_offset_ = 0;
	qreal TOPPASBase::z_value_ = 42.0;

  TOPPASBase::TOPPASBase(QWidget* parent):
      QMainWindow(parent),
      DefaultParamHandler("TOPPASBase"),
      clipboard_scene_(0)
  {
  	setWindowTitle("TOPPAS");
    setWindowIcon(QIcon(":/TOPPAS.png"));
		
    //prevents errors caused by too small width,height values
    setMinimumSize(400,400);

    // center main window
    setGeometry(
      (int)(0.1 * QApplication::desktop()->width()),
      (int)(0.1 * QApplication::desktop()->height()),
      (int)(0.8 * QApplication::desktop()->width()),
      (int)(0.8 * QApplication::desktop()->height())
      );

    // create dummy widget (to be able to have a layout), Tab bar and workspace
    QWidget* dummy = new QWidget(this);
    setCentralWidget(dummy);
    QVBoxLayout* box_layout = new QVBoxLayout(dummy);
    tab_bar_ = new TOPPASTabBar(dummy);
    tab_bar_->setWhatsThis("Tab bar<BR><BR>Close tabs through the context menu or by double-clicking them.");
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
    file->addAction("&New",this,SLOT(newPipeline()), Qt::CTRL+Qt::Key_N);
		file->addAction("Open &example file",this,SLOT(openExampleDialog()));
		file->addAction("&Open",this,SLOT(openFileDialog()), Qt::CTRL+Qt::Key_O);
    file->addAction("&Include",this,SLOT(includePipeline()), Qt::CTRL+Qt::Key_I);
    file->addAction("&Save",this,SLOT(savePipeline()), Qt::CTRL+Qt::Key_S);
    file->addAction("Save &As",this,SLOT(saveCurrentPipelineAs()), Qt::CTRL+Qt::SHIFT+Qt::Key_S);
		file->addAction("Refresh &parameters",this,SLOT(refreshParameters()), Qt::CTRL+Qt::SHIFT+Qt::Key_P);
    file->addAction("&Close",this,SLOT(closeFile()), Qt::CTRL+Qt::Key_W);
		file->addSeparator();
    file->addAction("&Load TOPPAS resource file",this,SLOT(loadPipelineResourceFile()));
    file->addAction("Save TOPPAS &resource file",this,SLOT(savePipelineResourceFile()));
    file->addSeparator();
    file->addAction("&Quit",qApp,SLOT(quit()));

    //Advanced menu
    //QMenu* advanced = new QMenu("&Advanced",this);
    //menuBar()->addMenu(advanced);
    //advanced->addAction("&Refresh definitions",this,SLOT(refreshDefinitions()), Qt::CTRL+Qt::Key_R);

		//Pipeline menu
		QMenu* pipeline = new QMenu("&Pipeline", this);
		menuBar()->addMenu(pipeline);
		pipeline->addAction("&Run (F5)",this,SLOT(runPipeline()));
		pipeline->addAction("&Abort",this,SLOT(abortPipeline()));

		//Windows menu
    QMenu * windows = new QMenu("&Windows", this);
    menuBar()->addMenu(windows);

		//Help menu
		QMenu* help = new QMenu("&Help", this);
		menuBar()->addMenu(help);
		QAction* action = help->addAction("OpenMS website",this,SLOT(showURL()));
		action->setData("http://www.OpenMS.de");
		action = help->addAction("TOPPAS tutorial (online)",this,SLOT(showURL()), Qt::Key_F1);
		action->setData("http://www-bs2.informatik.uni-tuebingen.de/services/OpenMS-release/html/TOPP_tutorial.html");
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

		//################## DEFAULTS #################
    //general
   	defaults_.setValue("preferences:default_path", ".", "Default path for loading and storing files.");
    defaults_.setValue("preferences:default_path_current", "true", "If the current path is preferred over the default path.");
		defaults_.setValidStrings("preferences:default_path_current",StringList::create("true,false"));
		defaults_.setValue("preferences:version","none","OpenMS version, used to check if the TOPPAS.ini is up-to-date");
				
  	defaultsToParam_();

  	//load param file
    loadPreferences();

		//################## Dock widgets #################
    //TOPP tools window
    QDockWidget* topp_tools_bar = new QDockWidget("TOPP", this);
    addDockWidget(Qt::LeftDockWidgetArea, topp_tools_bar);
    tools_tree_view_ = createTOPPToolsTreeWidget(topp_tools_bar);
    topp_tools_bar->setWidget(tools_tree_view_);    
		connect (tools_tree_view_, SIGNAL(itemDoubleClicked(QTreeWidgetItem*,int)), this, SLOT(insertNewVertexInCenter_(QTreeWidgetItem*)));
		windows->addAction(topp_tools_bar->toggleViewAction());

		//log window
		QDockWidget* log_bar = new QDockWidget("Log", this);
		addDockWidget(Qt::BottomDockWidgetArea, log_bar);
		log_ = new QTextEdit(log_bar);
		log_->setReadOnly(true);
		log_bar->setWidget(log_);
		log_bar->hide();
    //windows->addAction("&Show log window",log_bar,SLOT(show()));
		windows->addAction(log_bar->toggleViewAction());
		
		//set current path
		current_path_ = param_.getValue("preferences:default_path");
		
		//set temporary path
		tmp_path_ = QDir::tempPath() + QDir::separator();
		
  	//update the menu
  	updateMenu();
	}

  TOPPASBase::~TOPPASBase()
  {
  	savePreferences();
  }

  //static
  TOPPASTreeView* TOPPASBase::createTOPPToolsTreeWidget(QWidget* parent_widget)
  {
    TOPPASTreeView* tools_tree_view = new TOPPASTreeView(parent_widget);
    tools_tree_view->setWhatsThis("TOPP tools list<BR><BR>All available TOPP tools are shown here.");
    tools_tree_view->setColumnCount(1);
    QStringList header_labels;
    header_labels.append(QString("TOPP tools"));
    tools_tree_view->setHeaderLabels(header_labels);

    QTreeWidgetItem* item = new QTreeWidgetItem((QTreeWidget*)0);
    item->setText(0, "<Input files>");
    tools_tree_view->addTopLevelItem(item);
    item = new QTreeWidgetItem((QTreeWidget*)0);
    item->setText(0, "<Output files>");
    tools_tree_view->addTopLevelItem(item);
    item = new QTreeWidgetItem((QTreeWidget*)0);
    item->setText(0, "<Merger>");
    tools_tree_view->addTopLevelItem(item);

    //Param category_param = param_.copy("tool_categories:", true);

    ToolListType tools_list = ToolHandler::getTOPPToolList(true);
    ToolListType util_list = ToolHandler::getUtilList();
    // append utils
    for (ToolListType::Iterator it = util_list.begin(); it != util_list.end(); ++it)
    {
      it->second.category = "Utils";
      tools_list.insert(*it);
    }

    // any tool without a category gets into "unassigned" bin
    for (ToolListType::Iterator it = tools_list.begin(); it != tools_list.end(); ++it)
    {
      if (it->second.category.trim() =="") it->second.category = "Unassigned";
    }

    QSet<QString> category_set;
    for (ToolListType::ConstIterator it = tools_list.begin(); it != tools_list.end(); ++it)
    {
      category_set << String(it->second.category).toQString();
    }
    QStringList category_list = category_set.toList();
    qSort(category_list);
    Map<QString,QTreeWidgetItem*> category_map;

    foreach (const QString& category, category_list)
    {
      item = new QTreeWidgetItem((QTreeWidget*)0);
      item->setText(0, category);
      tools_tree_view->addTopLevelItem(item);
      category_map[category] = item;
    }

    QTreeWidgetItem* parent_item;
    for (ToolListType::iterator it = tools_list.begin(); it != tools_list.end(); ++it)
    {
      item = new QTreeWidgetItem(category_map[it->second.category.toQString()]);
      item->setText(0, it->first.toQString());
      parent_item = item;
      StringList types = ToolHandler::getTypes(it->first);
      for (StringList::iterator types_it = types.begin(); types_it != types.end(); ++types_it)
      {
        item = new QTreeWidgetItem(parent_item);
        item->setText(0, types_it->toQString());
      }
    }
    tools_tree_view->resizeColumnToContents(0);
    return tools_tree_view;
  }

	void TOPPASBase::loadFiles(const StringList& list, QSplashScreen* splash_screen)
  {
    for (StringList::const_iterator it=list.begin(); it!=list.end(); ++it)
    {
			splash_screen->showMessage((String("Loading file: ") + *it).toQString());
			splash_screen->repaint();
			QApplication::processEvents();
      addTOPPASFile(*it);
    }
  }

	void TOPPASBase::refreshDefinitions()
	{
	
	}
	
	void TOPPASBase::openExampleDialog()
	{
		QString file_name = QFileDialog::getOpenFileName(this, tr("Open example workflow"),
													File::getOpenMSDataPath().toQString()
													+QDir::separator()+"examples"+QDir::separator()
													+"TOPPAS"+QDir::separator(),
													tr("TOPPAS pipelines (*.toppas)"));
		
    addTOPPASFile(file_name);
	}
	
	void TOPPASBase::openFileDialog()
  {
		QString file_name = QFileDialog::getOpenFileName(this, tr("Open workflow"), current_path_.toQString(), tr("TOPPAS pipelines (*.toppas)"));
		
    addTOPPASFile(file_name);
  }
 	
  void TOPPASBase::includePipeline()
	{
		QString file_name = QFileDialog::getOpenFileName(this, tr("Include workflow"), current_path_.toQString(), tr("TOPPAS pipelines (*.toppas)"));		
    addTOPPASFile(file_name, false);
	}

  void TOPPASBase::addTOPPASFile(const String& file_name, bool in_new_window)
	{
		if (file_name != "")
		{
      if (!file_name.toQString().endsWith(".toppas", Qt::CaseInsensitive))
			{
				std::cerr << "The file '" << file_name << "' is not a .toppas file" << std::endl;
				return;
			}
		
			TOPPASScene* scene = 0;
			if (in_new_window)
			{
				TOPPASWidget* tw = new TOPPASWidget(Param(), ws_, tmp_path_);
				showAsWindow_(tw, File::basename(file_name));
				scene = tw->getScene();
				scene->load(file_name);
        connect(scene, SIGNAL(saveMe()), this, SLOT(savePipeline()));
        connect(scene, SIGNAL(selectionCopied(TOPPASScene*)), this, SLOT(saveToClipboard(TOPPASScene*)));
        connect(scene, SIGNAL(requestClipboardContent()), this, SLOT(sendClipboardContent()));
        connect(scene, SIGNAL(mainWindowNeedsUpdate()), this, SLOT(updateMenu()));
        connect(scene, SIGNAL(openInTOPPView(QVector<QStringList>)), this, SLOT(openFilesInTOPPView(QVector<QStringList>)));
      }
			else
			{
				if (!activeWindow_())
				{
					return;
				}
				TOPPASScene* tmp_scene = new TOPPASScene(0, QDir::tempPath()+QDir::separator(), false);
				tmp_scene->load(file_name);
				scene = activeWindow_()->getScene();
				scene->include(tmp_scene);
				delete tmp_scene;
			}
			
			//connect signals/slots for log messages
			for (TOPPASScene::VertexIterator it = scene->verticesBegin(); it != scene->verticesEnd(); ++it)
			{
				TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(*it);
				if (tv)
				{
					connect (tv, SIGNAL(toolStarted()), this, SLOT(toolStarted()));
					connect (tv, SIGNAL(toolFinished()), this, SLOT(toolFinished()));
					connect (tv, SIGNAL(toolCrashed()), this, SLOT(toolCrashed()));
					connect (tv, SIGNAL(toolFailed()), this, SLOT(toolFailed()));
					connect (tv, SIGNAL(toppOutputReady(const QString&)), this, SLOT(updateTOPPOutputLog(const QString&)));
					continue;
				}
				TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(*it);
				if (oflv)
				{
					connect (oflv, SIGNAL(outputFileWritten(const String&)), this, SLOT(outputVertexFinished(const String&)));
					continue;
				}
			}
		}
	}

  void TOPPASBase::newPipeline()
  {
  	TOPPASWidget* tw = new TOPPASWidget(Param(), ws_, tmp_path_);
		TOPPASScene* ts = tw->getScene();
		connect (ts, SIGNAL(selectionCopied(TOPPASScene*)), this, SLOT(saveToClipboard(TOPPASScene*)));
		connect (ts, SIGNAL(requestClipboardContent()), this, SLOT(sendClipboardContent()));
    connect (ts, SIGNAL(saveMe()), this, SLOT(savePipeline()));
		connect (ts, SIGNAL(mainWindowNeedsUpdate()), this, SLOT(updateMenu()));
  	showAsWindow_(tw, "(Untitled)");
  }
	
  void TOPPASBase::savePipeline()
	{
		TOPPASWidget* w = 0;
		QObject* sendr = QObject::sender();
		QAction* save_button_clicked = qobject_cast<QAction*>(sendr);
		
		if (!save_button_clicked)
		{
			// scene has requested to be saved
			TOPPASScene* ts = qobject_cast<TOPPASScene*>(sendr);
			if (ts && ts->views().size() > 0)
			{
				w = qobject_cast<TOPPASWidget*>(ts->views().first());
			}
		}
		else
		{
			w = activeWindow_();
		}
		
		if (!w)
		{
			return;
		}
		
		QString file_name = w->getScene()->getSaveFileName().toQString();
		if (file_name != "")
		{
      // accept also upper case TOPPAS extensions, since
      // we also support them while loading
      if (!file_name.endsWith(".toppas", Qt::CaseInsensitive))
			{
				file_name += ".toppas";
			}
			w->getScene()->store(file_name);
		}
		else
		{
      QString savedFileName = TOPPASBase::savePipelineAs(w, current_path_.toQString());
      // update tab title
      if(savedFileName != "")
      {
        QString caption = File::basename(savedFileName).toQString();
        tab_bar_->setTabText(tab_bar_->currentIndex(), caption);
      }
		}
	}
	
  void TOPPASBase::saveCurrentPipelineAs()
  {
    TOPPASWidget* w = activeWindow_();
    QString file_name = TOPPASBase::savePipelineAs(w, current_path_.toQString());
    if (file_name != "")
    {
      QString caption = File::basename(file_name).toQString();
      tab_bar_->setTabText(tab_bar_->currentIndex(), caption);
    }
  }

  // static
  QString TOPPASBase::savePipelineAs(TOPPASWidget* w, QString current_path)
	{
		if (!w)
		{
      return "";
		}
		
    QString file_name = QFileDialog::getSaveFileName(w, tr("Save workflow"), current_path, tr("TOPPAS pipelines (*.toppas)"));
		if (file_name != "")
		{
      if (!file_name.endsWith(".toppas", Qt::CaseInsensitive))
			{
				file_name += ".toppas";
			}
			w->getScene()->store(file_name);
      QString caption = File::basename(file_name).toQString();
      w->setWindowTitle(caption);
		}
    return file_name;
	}
	
  void TOPPASBase::loadPipelineResourceFile()
	{
		TOPPASWidget* w = activeWindow_();
    TOPPASBase::loadPipelineResourceFile(w, current_path_.toQString());
	}

  // static
  QString TOPPASBase::loadPipelineResourceFile(TOPPASWidget* w, QString current_path)
  {
    if (!w)
    {
      return "";
    }
    TOPPASScene* scene = w->getScene();
    QString file_name = QFileDialog::getOpenFileName(w, tr("Load resource file"), current_path, tr("TOPPAS resource files (*.trf)"));
    if (file_name == "")
    {
      return "";
    }
    TOPPASResources resources;
    resources.load(file_name);
    scene->loadResources(resources);
    return file_name;
  }
	
  void TOPPASBase::savePipelineResourceFile()
	{
		TOPPASWidget* w = activeWindow_();
    TOPPASBase::savePipelineResourceFile(w, current_path_.toQString());
	}

  // static
  QString TOPPASBase::savePipelineResourceFile(TOPPASWidget* w, QString current_path)
  {
    if (!w)
    {
      return "";
    }
    TOPPASScene* scene = w->getScene();
    QString file_name = QFileDialog::getSaveFileName(w, tr("Save resource file"), current_path, tr("TOPPAS resource files (*.trf)"));
    if (file_name == "")
    {
      return "";
    }
    if (!file_name.endsWith(".trf"))
    {
      file_name += ".trf";
    }
    TOPPASResources resources;
    scene->createResources(resources);
    resources.store(file_name);
    return file_name;
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
    connect(tw,SIGNAL(pipelineDroppedOnWidget(const String&, bool)),this,SLOT(addTOPPASFile(const String&, bool)));
	  tw->setWindowTitle(caption.toQString());

		//add tab with id
  	static int window_counter = 1337;
    tw->setWindowId(window_counter);
    window_counter++;

    tab_bar_->addTab(caption.toQString(), tw->getWindowId());

    //connect slots and sigals for removing the widget from the bar, when it is closed
    //- through the menu entry
    //- through the tab bar
    //- through the MDI close button
    connect(tw,SIGNAL(aboutToBeDestroyed(int)),tab_bar_,SLOT(removeId(int)));

    tab_bar_->setCurrentId(tw->getWindowId());

		//show first window maximized (only visible windows are in the list)
		if (ws_->windowList().count()==0)
		{
			tw->showMaximized();
		}
		else
		{
			tw->show();
		}
		TOPPASScene* ts = tw->getScene();
    connect (ts, SIGNAL(entirePipelineFinished()), this, SLOT(showPipelineFinishedLogMessage()));
		connect (ts, SIGNAL(entirePipelineFinished()), this, SLOT(updateMenu()));
		connect (ts, SIGNAL(pipelineExecutionFailed()), this, SLOT(updateMenu()));
		ts->setSceneRect((tw->mapToScene(tw->rect())).boundingRect());
  }

  void TOPPASBase::closeEvent(QCloseEvent* event)
  {
		bool close = true;
		QList<QWidget*> all_windows = ws_->windowList();
		foreach (QWidget* w, all_windows)
		{
			bool close_this = qobject_cast<TOPPASWidget*>(w)->getScene()->saveIfChanged();
			if (!close_this)
			{
				close = false;
				break;
			}
		}
		if (close)
		{
			//ws_->closeAllWindows();
			event->accept();
		}
		else
		{
			event->ignore();
		}
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
      if (window->getWindowId() == id)
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
  		// if the window refused to close, do not remove tab
			if (!window->isVisible())
			{
				tab_bar_->removeId(id);
			}
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
    QApplication::processEvents();
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
      Int window_id = qobject_cast<TOPPASWidget*>(w)->getWindowId();
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
      try
      { // the file might be corrupt
    	  tmp.load(filename);
      }
      catch (...)
      {
        error = true;
      }

    	//apply preferences if they are of the current TOPPAS version
    	if(!error && tmp.exists("preferences:version") && tmp.getValue("preferences:version").toString()==VersionInfo::getVersion())
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
  			//reset parameters (they will be stored again when TOPPAS quits)
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

    Param save_param = param_.copy("preferences:");
    
    try
    {
      // TODO: if closing multiple TOPPAS instances simultaneously, we might write to this file concurrently
      //       thus destroying its integrity. Think about using boost filelocks
      //       see OpenMS/METADATA/DocumentIDTagger.h for example
      //       and also implement in TOPPView (and other GUI's which write to user directory)
      save_param.store(string(param_.getValue("PreferencesFile")));
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

		QLabel* label = new QLabel(dlg);
		label->setPixmap(QPixmap(":/TOPP_about.png"));
		grid->addWidget(label,0,0);

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
		QLabel* text_label = new QLabel(text,dlg);
		grid->addWidget(text_label,0,1,Qt::AlignTop | Qt::AlignLeft);

		//execute
		dlg->exec();
	}
	
  void TOPPASBase::updateMenu()
  {
  	TOPPASWidget* tw = activeWindow_();
  	TOPPASScene* ts = 0;
  	if (tw)
  	{
  		ts = tw->getScene();
  	}
  	
  	QList<QAction*> actions = this->findChildren<QAction*>("");
		for (int i=0; i<actions.count(); ++i)
		{
			QString text = actions[i]->text();
			
			if (text=="&Run (F5)")
			{
				bool show = false;
				if (ts && !(ts->isPipelineRunning()))
				{
					show = true;
				}
				actions[i]->setEnabled(show);
			}
			else if (text=="&Abort")
			{
				bool show = false;
				if (ts && ts->isPipelineRunning())
				{
					show = true;
				}
				actions[i]->setEnabled(show);
			}
			else if (text=="&Include")
			{
				bool show = ts;
				actions[i]->setEnabled(show);
			}
			else if (text=="&Load resource file")
			{
				bool show = ts;
				actions[i]->setEnabled(show);
			}
			else if (text=="Save &resource file")
			{
				bool show = ts;
				actions[i]->setEnabled(show);
			}
			else if (text=="&Save")
			{
				bool show = ts && ts->wasChanged();
				actions[i]->setEnabled(show);
			}
			else if (text=="Refresh &parameters")
			{
				bool show = ts && !(ts->isPipelineRunning());
				actions[i]->setEnabled(show);
			}
		}
		
		if (ts)
		{
			QString title = tw->windowTitle();
			bool asterisk_shown = title.startsWith("*");
			bool changed = ts->wasChanged();
			if (asterisk_shown ^ changed)
			{
				title = asterisk_shown ? title.right(title.size()-1) : QString("*") + title;
				tw->setWindowTitle(title);
				tab_bar_->setTabText(tab_bar_->currentIndex(), title);
			}
		}
  }

  void TOPPASBase::showLogMessage_(TOPPASBase::LogState state, const String& heading, const String& body)
  {
		//Compose current time string
		DateTime d = DateTime::now();

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
		log_->moveCursor(QTextCursor::End);
  }

	void TOPPASBase::keyPressEvent(QKeyEvent* e)
	{
 		if (e->key() == Qt::Key_F5)
		{
			TOPPASWidget* tw = activeWindow_();
			if (!tw)
			{	
				e->ignore();
				return;
			}
			TOPPASScene* ts = tw->getScene();
			ts->runPipeline();
			e->accept();
		}
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
	
	void TOPPASBase::insertNewVertex_(double x, double y, QTreeWidgetItem* item)
	{
		if (!activeWindow_() || !activeWindow_()->getScene() || !tools_tree_view_)
		{
			return;
		}
		
		TOPPASScene* scene = activeWindow_()->getScene();
		QTreeWidgetItem* current_tool = item ? item : tools_tree_view_->currentItem();
		String tool_name = String(current_tool->text(0));
		TOPPASVertex* tv = 0;
		
		if (tool_name == "<Input files>")
		{
			tv = new TOPPASInputFileListVertex();
		}
		else if (tool_name == "<Output files>")
		{
			tv = new TOPPASOutputFileListVertex();
			TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(tv);
			connect (oflv, SIGNAL(outputFileWritten(const String&)), this, SLOT(outputVertexFinished(const String&)));
			scene->connectOutputVertexSignals(oflv);
		}
		else if (tool_name == "<Merger>")
		{
			tv = new TOPPASMergerVertex();
		}
		else // node is a TOPP tool
		{	
			if (current_tool->childCount() > 0)
			{
				// category or tool name with types is selected (instead of a concrete type)
				return;
			}
			String tool_type;
			if (current_tool->parent() != 0 && current_tool->parent()->parent() != 0)
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
			
			tv = new TOPPASToolVertex(tool_name, tool_type, tmp_path_);
			TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(tv);
			connect (ttv, SIGNAL(toolStarted()), this, SLOT(toolStarted()));
			connect (ttv, SIGNAL(toolFinished()), this, SLOT(toolFinished()));
			connect (ttv, SIGNAL(toolCrashed()), this, SLOT(toolCrashed()));
			connect (ttv, SIGNAL(toolFailed()), this, SLOT(toolFailed()));
			connect (ttv, SIGNAL(toppOutputReady(const QString&)), this, SLOT(updateTOPPOutputLog(const QString&)));
			
			scene->connectToolVertexSignals(ttv);
		}
		
		scene->connectVertexSignals(tv);
		scene->addVertex(tv);
		tv->setPos(x,y);
		tv->setZValue(z_value_);
		z_value_ += 0.000001;
		scene->topoSort();
		scene->setChanged(true);
	}

	void TOPPASBase::runPipeline()
	{
		TOPPASWidget* w = activeWindow_();
		if (w)
		{
			w->getScene()->runPipeline();
		}
	}
	
	void TOPPASBase::abortPipeline()
	{
		TOPPASWidget* w = activeWindow_();
		if (w)
		{
			w->getScene()->abortPipeline();
		}
		updateMenu();
	}
	
	void TOPPASBase::toolStarted()
	{
		TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
		if (tv)
		{
			String text = tv->getName();
			String type = tv->getType();
			if (type != "")
			{
				text += " ("+type+")";
			}
			text += " started. Processing ...";
			
			showLogMessage_(LS_NOTICE, text, "");
		}
		updateMenu();
	}
	
	void TOPPASBase::toolFinished()
	{
		TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
		if (tv)
		{
			String text = tv->getName();
			String type = tv->getType();
			if (type != "")
			{
				text += " ("+type+")";
			}
			text += " finished!";
			
			showLogMessage_(LS_NOTICE, text, "");
		}
		updateMenu();
	}
	
	void TOPPASBase::toolCrashed()
	{
		TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
		if (tv)
		{
			String text = tv->getName();
			String type = tv->getType();
			if (type != "")
			{
				text += " ("+type+")";
			}
			text += " crashed!";
			
			showLogMessage_(LS_ERROR, text, "");
		}
		updateMenu();
	}
	
	void TOPPASBase::toolFailed()
	{
		TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
		if (tv)
		{
			String text = tv->getName();
			String type = tv->getType();
			if (type != "")
			{
				text += " ("+type+")";
			}
			text += " failed!";
			
			showLogMessage_(LS_ERROR, text, "");
		}
		updateMenu();
	}
	
	void TOPPASBase::outputVertexFinished(const String& file)
	{
		String text = "Output file '"+file+"' written.";
		showLogMessage_(LS_NOTICE, text, "");
	}
	
	void TOPPASBase::updateTOPPOutputLog(const QString& out)
	{
		TOPPASToolVertex* sender = qobject_cast<TOPPASToolVertex*>(QObject::sender());
		if (!sender)
		{
			return;
		}
    /*
		QString text = (sender->getName()).toQString();
		if (sender->getType() != "")
		{
			text += " ("+(sender->getType()).toQString()+")";
		}
		text += ":\n" + out;
    */
    QString text = out; // shortened version for now (if we reintroduce simultaneous tool execution, 
                        // we need to rethink this (probably only trigger this slot when tool 100% finished)

		
		//show log if there is output
		qobject_cast<QWidget*>(log_->parent())->show();

		//update log_
		log_->append(text);
	}
	
  void TOPPASBase::showPipelineFinishedLogMessage()
	{
		showLogMessage_(LS_NOTICE, "Entire pipeline execution finished!", "");
	}
	
	void TOPPASBase::insertNewVertexInCenter_(QTreeWidgetItem* item)
	{
		if (!activeWindow_() || !activeWindow_()->getScene() || !tools_tree_view_ || !tools_tree_view_->currentItem())
		{
			return;
		}
		
		QPointF insert_pos = activeWindow_()->mapToScene(QPoint((activeWindow_()->width()/2.0)+(qreal)(5*node_offset_),(activeWindow_()->height()/2.0)+(qreal)(5*node_offset_)));
		insertNewVertex_(insert_pos.x(), insert_pos.y(), item);
		node_offset_ = (node_offset_+1) % 10;
	}
	
	void TOPPASBase::saveToClipboard(TOPPASScene* scene)
	{
    if (clipboard_scene_ != 0)
		{
      delete clipboard_scene_;
      clipboard_scene_ = 0;
		}
    clipboard_scene_ = scene;
	}

	void TOPPASBase::sendClipboardContent()
	{
		TOPPASScene* sndr = qobject_cast<TOPPASScene*>(QObject::sender());
		if (sndr != 0)
		{
      sndr->setClipboard(clipboard_scene_);
		}
	}
	
  void TOPPASBase::refreshParameters()
	{
    TOPPASWidget* w = activeWindow_();
    QString file_name = TOPPASBase::refreshPipelineParameters(w, current_path_.toQString());
    if (file_name != "")
    {
      QString caption = File::basename(file_name).toQString();
      tab_bar_->setTabText(tab_bar_->currentIndex(), caption);
    }
	}

  // static
  QString TOPPASBase::refreshPipelineParameters(TOPPASWidget* tw, QString current_path)
  {
    TOPPASScene* ts = 0;
    if (tw)
    {
      ts = tw->getScene();
    }
    if (!ts)
    {
      return "";
    }

    if (!ts->refreshParameters())
    {
      QMessageBox::information(tw, tr("Nothing to be done"),
                               tr("The parameters of the tools used in this workflow have not changed."));
      return "";
    }

    ts->setChanged(true);
    int ret = QMessageBox::information(tw, "Parameters updated!",
                                       "The parameters of some tools in this workflow have changed. Do you want to save these changes now?",
                                       QMessageBox::Save | QMessageBox::Cancel);
    if (ret == QMessageBox::Save)
    {
      QString file_name = TOPPASBase::savePipelineAs(tw, current_path);
      return file_name;
    }

    return "";
  }

  void TOPPASBase::openFilesInTOPPView(QVector<QStringList> all_files)
  {    
    foreach (const QStringList& files, all_files)
    {
      if (files.size() > 0)
      {
        QProcess* p = new QProcess();
        p->setProcessChannelMode(QProcess::ForwardedChannels);
        QString toppview_executable;
        toppview_executable = "TOPPView";        
        QStringList arg = files;

        if (files.size() > 1)
        {
          // ask user how to open multiple files
          if ( !QMessageBox::question(
                  this,
                  tr("Open in separate windows? -- TOPPAS"),
                  tr("How do you want to open the output files?"),
                  tr("&Single window"), tr("&Separate windows"),
                  QString::null, 0, 1 ) )
          {
            arg = files.join(" + ").split(" ", QString::SkipEmptyParts);
          }
        }

        p->start(toppview_executable, arg);
        if (!p->waitForStarted())
        {
          // execution failed
          std::cerr << p->errorString().toStdString() << std::endl;
#if defined(Q_WS_MAC)
          std::cerr << "Please check if TOPPAS and TOPPView are located in the same directory" << std::endl;
#endif
        }
      }
    }
  }

} //namespace OpenMS


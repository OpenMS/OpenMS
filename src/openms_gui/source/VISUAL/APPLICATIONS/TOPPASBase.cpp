// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/APPLICATIONS/TOPPASBase.h>

#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/VISUAL/EnhancedWorkspace.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/LogWindow.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASResources.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASSplitterVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASWidget.h>

#include <map>

//Qt
#include <QApplication>
#include <QCloseEvent>
#include <QDesktopServices>
#include <QNetworkAccessManager>
#include <QNetworkProxy>
#include <QNetworkProxyFactory>
#include <QNetworkReply>
#include <QSvgGenerator>
#include <QTextStream>
#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QMap>
#include <QtCore/QSet>
#include <QtCore/QSettings>
#include <QtCore/QUrl>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QInputDialog>
#include <QtWidgets/QLabel>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QListWidgetItem>
#include <QtWidgets/QMdiSubWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSplashScreen>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QToolTip>
#include <QtWidgets/QTreeWidget>
#include <QtWidgets/QTreeWidgetItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWhatsThis>
#include <utility>
#include <OpenMS/VISUAL/TOPPASOutputFolderVertex.h>


using namespace std;
using namespace OpenMS;

namespace OpenMS
{
  using namespace Internal;

  int TOPPASBase::node_offset_ = 0;
  qreal TOPPASBase::z_value_ = 42.0;

  TOPPASBase::TOPPASBase(QWidget* parent) :
    QMainWindow(parent),
    DefaultParamHandler("TOPPASBase"),
    clipboard_scene_(nullptr)
  {

    setWindowTitle("TOPPAS");
    setWindowIcon(QIcon(":/TOPPAS.png"));

    //prevents errors caused by too small width,height values
    setMinimumSize(400, 400);

    // center main window
    setGeometry(
      (int)(0.1 * QGuiApplication::primaryScreen()->geometry().width()),
      (int)(0.1 * QGuiApplication::primaryScreen()->geometry().height()),
      (int)(0.8 * QGuiApplication::primaryScreen()->geometry().width()),
      (int)(0.8 * QGuiApplication::primaryScreen()->geometry().height())
    );

    // create dummy widget (to be able to have a layout), Tab bar and workspace
    QWidget* dummy = new QWidget(this);
    setCentralWidget(dummy);
    QVBoxLayout* box_layout = new QVBoxLayout(dummy);
    tab_bar_ = new EnhancedTabBar(dummy);
    tab_bar_->setWhatsThis("Tab bar<BR><BR>Close tabs through the context menu or by double-clicking them.");
    tab_bar_->addTab("dummy", 1336);
    tab_bar_->setMinimumSize(tab_bar_->sizeHint());
    tab_bar_->removeId(1336);
    //connect slots and signals for selecting spectra
    connect(tab_bar_, &EnhancedTabBar::currentIdChanged, this, &TOPPASBase::focusByTab);
    connect(tab_bar_, &EnhancedTabBar::closeRequested, this, &TOPPASBase::closeByTab);

    box_layout->addWidget(tab_bar_);
    ws_ = new EnhancedWorkspace(dummy);
    connect(ws_, &EnhancedWorkspace::subWindowActivated, this, &TOPPASBase::updateTabBar);
    connect(ws_, &EnhancedWorkspace::subWindowActivated, this, &TOPPASBase::updateMenu);

    box_layout->addWidget(ws_);

    //################## MENUS #################
    // File menu
    QMenu* file = new QMenu("&File", this);
    menuBar()->addMenu(file);
    file->addAction("&New", this, SLOT(newPipeline()), Qt::CTRL | Qt::Key_N);
    file->addAction("&Open", this, SLOT(openFilesByDialog()), Qt::CTRL + Qt::Key_O);
    file->addAction("Open &example file", this, SLOT(openExampleDialog()), Qt::CTRL + Qt::Key_E);
    file->addAction("&Include", this, SLOT(includePipeline()), Qt::CTRL + Qt::Key_I);
    //file->addAction("Online &Repository", this, SLOT(openOnlinePipelineRepository()), Qt::CTRL + Qt::Key_R);
    file->addAction("&Save", this, SLOT(savePipeline()), Qt::CTRL + Qt::Key_S);
    file->addAction("Save &As", this, SLOT(saveCurrentPipelineAs()), Qt::CTRL | Qt::SHIFT | Qt::Key_S);
    file->addAction("E&xport as image", this, SLOT(exportAsImage()));
    file->addAction("Refresh &parameters", this, SLOT(refreshParameters()), Qt::CTRL | Qt::SHIFT | Qt::Key_P);
    file->addAction("&Close pipeline", this, SLOT(closeFile()), Qt::CTRL | Qt::Key_W);

    file->addSeparator();
    // Recent files
    file->addMenu(recent_files_menu_.getMenu()); // updates automatically via RecentFilesMenu class, since this is just a pointer
    connect(&recent_files_menu_, &RecentFilesMenu::recentFileClicked, [this](const String& filename) { addTOPPASFile(filename, true);});

    file->addSeparator();
    file->addAction("&Load TOPPAS resource file", this, SLOT(loadPipelineResourceFile()));
    file->addAction("Sa&ve TOPPAS resource file", this, SLOT(savePipelineResourceFile()));
    file->addSeparator();
    file->addAction("&Quit", qApp, SLOT(quit()));

    //Pipeline menu
    QMenu* pipeline = new QMenu("&Pipeline", this);
    menuBar()->addMenu(pipeline);
    pipeline->addAction("&Run (F5)", this, SLOT(runPipeline()));
    pipeline->addAction("&Abort", this, SLOT(abortPipeline()));

    //Windows menu
    QMenu* windows = new QMenu("&Windows", this);
    menuBar()->addMenu(windows);

    //Help menu
    QMenu* help = new QMenu("&Help", this);
    menuBar()->addMenu(help);
    QAction* action = help->addAction("OpenMS website", this, SLOT(showURL()));
    action->setData("http://www.OpenMS.de");
    action = help->addAction("TOPPAS tutorial", this, SLOT(showURL()), Qt::Key_F1);
    action->setData(String("html/TOPPAS_tutorial.html").toQString());

    help->addSeparator();
    help->addAction("&About", this, SLOT(showAboutDialog()));


    //create status bar
    message_label_ = new QLabel(statusBar());
    statusBar()->addWidget(message_label_, 1);

    //################## TOOLBARS #################
    //create toolbars and connect signals

    //--Basic tool bar--
    //tool_bar_ = addToolBar("Basic tool bar");
    //tool_bar_->show();

    //################## DEFAULTS #################
    //general
    defaults_.setValue("preferences:default_path", ".", "Default path for loading and storing files.");
    defaults_.setValue("preferences:default_path_current", "true", "If the current path is preferred over the default path.");
    defaults_.setValidStrings("preferences:default_path_current", {"true","false"});
    defaults_.setValue("preferences:version", "none", "OpenMS version, used to check if the TOPPAS.ini is up-to-date");
    subsections_.emplace_back("preferences:RecentFiles");
    defaultsToParam_();

    //load param file
    loadPreferences();

    //################## Dock widgets #################
    //TOPP tools window
    QDockWidget* topp_tools_bar = new QDockWidget("TOPP", this);
    topp_tools_bar->setObjectName("TOPP_tools_bar");
    addDockWidget(Qt::LeftDockWidgetArea, topp_tools_bar);
    QWidget* frame = new QWidget(topp_tools_bar);
    auto frame_layout = new QVBoxLayout(frame);
    //frame->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Maximum);
    tools_tree_view_ = createTOPPToolsTreeWidget();
    tools_filter_ = new QLineEdit();
    tools_expand_all_ = new QPushButton("expand all");
    tools_collapse_all_ = new QPushButton("collapse all");
    frame_layout->addWidget(new QLabel("Filter: "));
    frame_layout->addWidget(tools_filter_);
    frame_layout->addWidget(tools_expand_all_);
    frame_layout->addWidget(tools_collapse_all_);
    frame_layout->addWidget(tools_tree_view_);
    topp_tools_bar->setWidget(frame);
    connect(tools_expand_all_, &QPushButton::clicked, tools_tree_view_, &TOPPASTreeView::expandAll);
    connect(tools_collapse_all_, &QPushButton::clicked, tools_tree_view_, &TOPPASTreeView::collapseAll);
    connect(tools_tree_view_, &QTreeWidget::itemDoubleClicked, this, &TOPPASBase::insertNewVertexInCenter_);
    connect(tools_filter_, &QLineEdit::textChanged, this, &TOPPASBase::filterToolTree_);
    windows->addAction(topp_tools_bar->toggleViewAction());

    //log window
    QDockWidget* log_bar = new QDockWidget("Log", this);
    log_bar->setObjectName("log_bar");
    addDockWidget(Qt::BottomDockWidgetArea, log_bar);
    log_ = new LogWindow(log_bar);
    log_->setMaxLength(1e7); // limit to 10 mio characters, and trim to 5 mio upon reaching this limit
    log_bar->setWidget(log_);
    log_bar->hide();
    //windows->addAction("&Show log window",log_bar,SLOT(show()));
    windows->addAction(log_bar->toggleViewAction());

    //workflow description window
    QDockWidget* description_bar = new QDockWidget("Workflow Description", this);
    description_bar->setObjectName("workflow_description_bar");
    addDockWidget(Qt::RightDockWidgetArea, description_bar);
    desc_ = new QTextEdit(description_bar);
    desc_->setTextColor(Qt::black);
    desc_->setText("... put your workflow description here ...");
    desc_->setTextColor(Qt::black);
    desc_->document()->setDefaultFont(QFont("Arial", 12));
    description_bar->setWidget(desc_);
    windows->addAction(description_bar->toggleViewAction());
    connect(desc_, SIGNAL(textChanged()), this, SLOT(descriptionUpdated_()));

    // set current path
    current_path_ = param_.getValue("preferences:default_path").toString();

    // set & create temporary path -- make sure its a new subdirectory, as it will be deleted later
    QString new_tmp_dir = File::getUniqueName(false).toQString();
    QDir qd(File::getTempDirectory().toQString());
    qd.mkdir(new_tmp_dir);
    qd.cd(new_tmp_dir);
    tmp_path_ = qd.absolutePath();

/* 
     QT5 replace with QWebEngine
    // online browser
    webview_ = new QWebView(parent);
    webview_->page()->setLinkDelegationPolicy(QWebPage::DelegateAllLinks); // now linkClicked() is emitted
    connect((webview_->page()), SIGNAL(linkClicked(const QUrl &)), this, SLOT(downloadTOPPASfromHomepage_(const QUrl &)));
*/

    network_manager_ = new QNetworkAccessManager(this);
    connect(network_manager_, SIGNAL(finished(QNetworkReply*)), this, SLOT(toppasFileDownloaded_(QNetworkReply*)));

    // update the menu
    updateMenu(); 

    QSettings settings("OpenMS", "TOPPAS");
    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());
  }

  void TOPPASBase::filterToolTree_()
  {
    tools_tree_view_->filter(tools_filter_->text());
  }

  TOPPASBase::~TOPPASBase()
  {
    savePreferences();
    // delete temporary files (TODO: make this a user dialog and ask - for later resume)
    // safety measure: only delete if subdirectory of Temp path; we do not want to delete / or c:
    if (String(tmp_path_).substitute("\\", "/").hasPrefix(File::getTempDirectory().substitute("\\", "/") + "/"))
    {
      File::removeDirRecursively(tmp_path_);
    }
  }

  void TOPPASBase::descriptionUpdated_()
  {
    if (!activeSubWindow_() || !activeSubWindow_()->getScene())
    {
      return;
    }
    //std::cerr << "changed to '" << String(desc_->toHtml()) << "'\n";
    activeSubWindow_()->getScene()->setChanged(true);
    activeSubWindow_()->getScene()->setDescription(desc_->toHtml());
  }

  void TOPPASBase::toppasFileDownloaded_(QNetworkReply* /* r */)
  {
/* QT5
    r->deleteLater();
    if (r->error() != QNetworkReply::NoError)
    {
      log_->appendNewHeader(LogWindow::LogState::CRITICAL, "Download failed", "Error '" + r->errorString() + "' while downloading TOPPAS file: '" + r->url().toString() + "'");
      return;
    }

    QByteArray data = r->readAll();

    QString proposed_filename;
    if (r->url().hasQueryItem("file"))
    {
      proposed_filename = r->url().queryItemValue("file");
    }
    else
    {
      proposed_filename = "Workflow.toppas";
      OPENMS_LOG_WARN << "The URL format of downloads from the TOPPAS Online-Repository has changed. Please notify developers!";
    }
    QString filename = QFileDialog::getSaveFileName(this, "Where to save the TOPPAS file?", this->current_path_.toQString() + "/" + proposed_filename, tr("TOPPAS (*.toppas)"));

    // check if the user clicked cancel, to avoid saving .toppas somewhere
    if (String(filename).trim().empty())
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "Download succeeded, but saving aborted by user!", "");
      return;
    }

    if (!filename.endsWith(".toppas", Qt::CaseInsensitive))
    {
      filename += ".toppas";
    }

    QFile file(filename);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
      log_->appendNewHeader(LogWindow::LogState::NOTICE, "Download succeeded. Cannot save the file. Try again with another filename and/or location!", "");
      return;
    }

    QTextStream out(&file);
    out << data;
    file.close();

    this->addTOPPASFile(filename);
    log_->appendNewHeader(LogWindow::LogState::NOTICE, "File successfully saved to '" + filename + "'.", "");
*/
  }
  
  void TOPPASBase::TOPPASreadyRead()
  {
    QNetworkReply::NetworkError ne = network_reply_->error();
    qint64 ba = network_reply_->bytesAvailable();
    OPENMS_LOG_DEBUG << "Error code (QNetworkReply::NetworkError): " << ne << "  bytes available: " << ba << std::endl;
    return;
  }

  void TOPPASBase::downloadTOPPASfromHomepage_(const QUrl& url)
  {
    if (url.toString().endsWith(QString(".toppas"), Qt::CaseInsensitive))
    {
      network_reply_ = network_manager_->get(QNetworkRequest(url));

      // debug
      connect(network_reply_, SIGNAL(readyRead()), this, SLOT(TOPPASreadyRead()));
      connect(network_reply_, SIGNAL(error(QNetworkReply::NetworkError code)), this, SLOT(TOPPASreadyRead()));
      connect(network_reply_, SIGNAL(finished()), this, SLOT(TOPPASreadyRead()));
      connect(network_reply_, SIGNAL(metaDataChanged()), this, SLOT(TOPPASreadyRead()));
      connect(network_reply_, SIGNAL(sslErrors(const QList<QSslError> & errors)), this, SLOT(TOPPASreadyRead()));
      // .. end debug

      log_->appendNewHeader(LogWindow::LogState::NOTICE, "Downloading file '" + url.toString() + "'. You will be notified once the download finished.", "");
      // webview_->close(); QT5 replace with QWebEngine
    }
    else
    {
      QMessageBox::warning(this, tr("Error"), tr("You can only click '.toppas' files on this page. No navigation is allowed!\n"));
      /* 
      replace with QT5 webengine
      webview_->setFocus(); 
      webview_->activateWindow();
      */
    }
  }

  void TOPPASBase::openOnlinePipelineRepository()
  {
/* QT5
    QUrl url = QUrl("http://www.OpenMS.de/TOPPASWorkflows/");

    static bool proxy_settings_checked = false;
    if (!proxy_settings_checked) // do only once because may take several seconds on windows
    {
      QNetworkProxy proxy;
      QUrl tmp_proxy_url(QString(getenv("http_proxy")));
      QUrl tmp_PROXY_url(QString(getenv("HTTP_PROXY")));
      QUrl proxy_url = tmp_proxy_url.isValid() ? tmp_proxy_url : tmp_PROXY_url;
      if (proxy_url.isValid())
      {
        QString hostname = proxy_url.host();
        int port = proxy_url.port();
        QString username = proxy_url.userName();
        QString password = proxy_url.password();
        proxy = QNetworkProxy(QNetworkProxy::HttpProxy, hostname, port, username, password);
      }
      else
      {
        QList<QNetworkProxy> proxies = QNetworkProxyFactory::systemProxyForQuery(url);
        if (!proxies.empty())
        {
          proxy = proxies.first();
        }
      }
      QNetworkProxy::setApplicationProxy(proxy); //no effect if proxy == QNetworkProxy()
      proxy_settings_checked = true;
    }

    // show something immediately, so the user does not stare at white screen while the URL is fetched
    webview_->setHtml("loading content ... ");
    webview_->show();
    // ... load the page in background
    webview_->load(url);
*/
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

    auto add_list_item = [&tools_tree_view](const QString& node_name, const QString& tool_tip)
      {
      QTreeWidgetItem* item = new QTreeWidgetItem(tools_tree_view);
      item->setText(0, node_name);
      item->setToolTip(0, tool_tip);
      tools_tree_view->addTopLevelItem(item);
    };
    add_list_item("<Input files>", "One or multiple input files, such as mzML or FASTA files from your local hard drive");
    add_list_item("<Output files>", "Sink for one or more output files, which are produced by a TOPP tool and which you want to keep for later.");
    add_list_item("<Output folder>", "Some TOPP tools write their output to a folder. Usually a fixed set of files, whose names cannot be set explicitly.");
    add_list_item("<Merger>", "Concatenate files from multiple input edges to a list and forward that list.");
    add_list_item("<Collector>", "Collect each single file from \na single input edge (for every time it runs)\nand then foward this list to the next tool (which is only invoked once)");
    add_list_item("<Splitter>", "Opposite of a collector.");

    //Param category_param = param_.copy("tool_categories:", true);

    ToolListType tools_list = ToolHandler::getTOPPToolList(true);

    // any tool without a category gets into "unassigned" bin
    for (ToolListType::iterator it = tools_list.begin(); it != tools_list.end(); ++it)
    {
      if (it->second.category.trim().empty())
        it->second.category = "Unassigned";
    }

    QSet<QString> category_set;
    for (ToolListType::const_iterator it = tools_list.begin(); it != tools_list.end(); ++it)
    {
      category_set << String(it->second.category).toQString();
    }

    QStringList category_list = category_set.values();
    std::sort(category_list.begin(), category_list.end());

    std::map<QString, QTreeWidgetItem*> category_map;

    for (const QString &category : category_list)
    {
      auto item = new QTreeWidgetItem((QTreeWidget*)nullptr);
      item->setText(0, category);
      tools_tree_view->addTopLevelItem(item);
      category_map[category] = item;
    }

    for (const auto& tool : tools_list)
    {
      auto item = new QTreeWidgetItem(category_map[tool.second.category.toQString()]);
      item->setText(0, tool.first.toQString());
      QTreeWidgetItem* parent_item = item;
      StringList types = ToolHandler::getTypes(tool.first);
      for (const auto& type : types)
      {
        item = new QTreeWidgetItem(parent_item);
        item->setText(0, type.toQString());
      }
    }
    tools_tree_view->resizeColumnToContents(0);
    return tools_tree_view;
  }

  void TOPPASBase::loadFiles(const StringList& list, QSplashScreen* splash_screen)
  {
    for (StringList::const_iterator it = list.begin(); it != list.end(); ++it)
    {
      splash_screen->showMessage((String("Loading file: ") + *it).toQString());
      splash_screen->repaint();
      QApplication::processEvents();
      addTOPPASFile(*it);
    }
  }

  void TOPPASBase::openExampleDialog()
  {
    QString file_name = QFileDialog::getOpenFileName(this, tr("Open example workflow"),
                                                     File::getOpenMSDataPath().toQString()
                                                     + QDir::separator() + "examples" + QDir::separator()
                                                     + "TOPPAS" + QDir::separator(),
                                                     tr("TOPPAS pipelines (*.toppas)"));

    addTOPPASFile(file_name);
  }

  void TOPPASBase::openFilesByDialog()
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
    if (file_name.empty()) return;

    if (!file_name.toQString().endsWith(".toppas", Qt::CaseInsensitive))
    {
      OPENMS_LOG_ERROR << "The file '" << file_name << "' is not a .toppas file" << std::endl;
      return;
    }
    
    recent_files_menu_.add(file_name);

    TOPPASWidget* asw = activeSubWindow_();
    TOPPASScene* scene = nullptr;
    if (in_new_window)
    {
      if (asw)
      {
        TOPPASWidget* uninitialized_window = window_(asw->getFirstWindowID());
        if (uninitialized_window && !uninitialized_window->getScene()->wasChanged())
          closeByTab(asw->getFirstWindowID());
      }
      TOPPASWidget* tw = new TOPPASWidget(Param(), ws_, tmp_path_);
      scene = tw->getScene();
      scene->load(file_name); // first load WF, including description etc
      showAsWindow_(tw, File::basename(file_name)); // show it
    }
    else
    {
      if (!activeSubWindow_()) return;

      TOPPASScene* tmp_scene = new TOPPASScene(nullptr, this->tmp_path_.toQString(), false);
      tmp_scene->load(file_name);
      scene = activeSubWindow_()->getScene();
      scene->include(tmp_scene);
      delete tmp_scene;
    }

    //connect signals/slots for log messages
    for (TOPPASScene::VertexIterator it = scene->verticesBegin(); it != scene->verticesEnd(); ++it)
    {
      TOPPASToolVertex* tv = dynamic_cast<TOPPASToolVertex*>(*it);
      if (tv)
      {
        connect(tv, SIGNAL(toolStarted()), this, SLOT(toolStarted()));
        connect(tv, SIGNAL(toolFinished()), this, SLOT(toolFinished()));
        connect(tv, SIGNAL(toolCrashed()), this, SLOT(toolCrashed()));
        connect(tv, SIGNAL(toolFailed()), this, SLOT(toolFailed()));
        connect(tv, SIGNAL(toolFailed(const QString &)), this, SLOT(updateTOPPOutputLog(const QString &)));
        // already done in ToppasScene:
        //connect (tv, SIGNAL(toppOutputReady(const QString&)), this, SLOT(updateTOPPOutputLog(const QString&)));
        continue;
      }

      TOPPASMergerVertex* tmv = dynamic_cast<TOPPASMergerVertex*>(*it);
      if (tmv)
      {
        connect(tmv, SIGNAL(mergeFailed(const QString)), this, SLOT(updateTOPPOutputLog(const QString &)));
        continue;
      }

      TOPPASOutputFileListVertex* oflv = dynamic_cast<TOPPASOutputFileListVertex*>(*it);
      if (oflv)
      {
        connect(oflv, SIGNAL(outputFileWritten(const String &)), this, SLOT(outputVertexFinished(const String &)));
        continue;
      }
    }
  }

  void TOPPASBase::newPipeline()
  {
    TOPPASWidget* tw = new TOPPASWidget(Param(), ws_, tmp_path_);
    showAsWindow_(tw, "(Untitled)");
  }

  void TOPPASBase::savePipeline()
  {
    TOPPASWidget* w = nullptr;
    QObject* sendr = QObject::sender();
    QAction* save_button_clicked = dynamic_cast<QAction*>(sendr);

    if (!save_button_clicked)
    {
      // scene has requested to be saved
      TOPPASScene* ts = dynamic_cast<TOPPASScene*>(sendr);
      if (ts && !ts->views().empty())
      {
        w = dynamic_cast<TOPPASWidget*>(ts->views().first());
      }
    }
    else
    {
      w = activeSubWindow_();
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
      if (!w->getScene()->store(file_name))
      {
        QMessageBox::warning(this, tr("Error"),
                             tr("Unable to save current pipeline. Possible reason: Invalid edges due to parameter refresh."));
      }

    }
    else
    {
      QString savedFileName = TOPPASBase::savePipelineAs(w, current_path_.toQString());
      // update tab title
      if (savedFileName != "")
      {
        tab_bar_->setTabText(File::basename(savedFileName).toQString());
      }
    }
  }

  void TOPPASBase::saveCurrentPipelineAs()
  {
    TOPPASWidget* w = activeSubWindow_();
    QString file_name = TOPPASBase::savePipelineAs(w, current_path_.toQString());
    if (file_name != "")
    {
      tab_bar_->setTabText(File::basename(file_name).toQString());
    }
  }

  // static
  QString TOPPASBase::savePipelineAs(TOPPASWidget* w, const QString& current_path)
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
      if (!w->getScene()->store(file_name))
      {
        QMessageBox::warning(nullptr, tr("Error"),
                             tr("Unable to save current pipeline. Possible reason: Invalid edges due to parameter refresh."));
      }
      QString caption = File::basename(file_name).toQString();
      w->setWindowTitle(caption);
    }
    return file_name;
  }

  void TOPPASBase::exportAsImage()
  {
    TOPPASWidget* w = activeSubWindow_();
    TOPPASScene* s = w->getScene();

    QString cp = current_path_.toQString();
    QString file_name = QFileDialog::getSaveFileName(w, tr("Save image"), cp, tr("Images (*.svg *.png *.jpg)"));
    if (file_name == "")
    {
      return;
    }
    if (!file_name.endsWith(".svg", Qt::CaseInsensitive) &&
        !file_name.endsWith(".png", Qt::CaseInsensitive) &&
        !file_name.endsWith(".jpg", Qt::CaseInsensitive))
    {
      file_name += ".svg";
    }
    bool svg = file_name.endsWith(".svg");

    QRectF items_bounding_rect = s->itemsBoundingRect();
    qreal wh_proportion = (qreal)(items_bounding_rect.width()) / (qreal)(items_bounding_rect.height());
    bool w_larger_than_h = wh_proportion > 1;
    qreal x1 = 0;
    qreal y1 = 0;
    qreal x2, y2;

    qreal small_edge_length = svg ? 500 : 4000;

    if (w_larger_than_h)
    {
      x2 = wh_proportion * small_edge_length;
      y2 = small_edge_length;
    }
    else
    {
      x2 = small_edge_length;
      y2 = (1.0 / wh_proportion) * small_edge_length;
    }
    qreal width = x2 - x1;
    qreal height = y2 - y1;

    if (svg)
    {
      QSvgGenerator svg_gen;
      svg_gen.setFileName(file_name);
      svg_gen.setSize(QSize(width, height));
      svg_gen.setViewBox(QRect(x1, y1, x2, y2));
      svg_gen.setTitle(tr("Title (TBD)"));
      svg_gen.setDescription(tr("Description (TBD)"));
      QPainter painter(&svg_gen);
      s->render(&painter, QRectF(), items_bounding_rect);
    }
    else
    {
      QImage img(width, height, QImage::Format_RGB32);
      img.fill(QColor(Qt::white).rgb());
      QPainter painter(&img);
      s->render(&painter, QRectF(), items_bounding_rect);
      img.save(file_name);
    }
  }

  void TOPPASBase::loadPipelineResourceFile()
  {
    TOPPASWidget* w = activeSubWindow_();
    TOPPASBase::loadPipelineResourceFile(w, current_path_.toQString());
  }

  // static
  QString TOPPASBase::loadPipelineResourceFile(TOPPASWidget* w, const QString& current_path)
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
    TOPPASWidget* w = activeSubWindow_();
    TOPPASBase::savePipelineResourceFile(w, current_path_.toQString());
  }

  // static
  QString TOPPASBase::savePipelineResourceFile(TOPPASWidget* w, const QString& current_path)
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
    ws_->addSubWindow(tw);
    tw->showMaximized();
    connect(tw, SIGNAL(sendStatusMessage(std::string, OpenMS::UInt)), this, SLOT(showStatusMessage(std::string, OpenMS::UInt)));
    connect(tw, SIGNAL(sendCursorStatus(double, double)), this, SLOT(showCursorStatus(double, double)));
    connect(tw, SIGNAL(toolDroppedOnWidget(double, double)), this, SLOT(insertNewVertex_(double, double)));
    connect(tw, SIGNAL(pipelineDroppedOnWidget(const String &, bool)), this, SLOT(addTOPPASFile(const String &, bool)));
    tw->setWindowTitle(caption.toQString());

    tw->addToTabBar(tab_bar_, caption, true);

    //show first window maximized (only visible windows are in the list)
    if (ws_->subWindowList().count() == 0)
    {
      tw->showMaximized();
    }
    else
    {
      tw->show();
    }
    TOPPASScene* scene = tw->getScene();
    connect(scene, SIGNAL(saveMe()), this, SLOT(savePipeline()));
    connect(scene, SIGNAL(selectionCopied(TOPPASScene*)), this, SLOT(saveToClipboard(TOPPASScene*)));
    connect(scene, SIGNAL(requestClipboardContent()), this, SLOT(sendClipboardContent()));
    connect(scene, SIGNAL(mainWindowNeedsUpdate()), this, SLOT(updateMenu()));
    connect(scene, SIGNAL(openInTOPPView(QStringList)), this, SLOT(openFilesInTOPPView(QStringList)));
    connect(scene, SIGNAL(messageReady(const QString &)), this, SLOT(updateTOPPOutputLog(const QString &)));
    connect(scene, SIGNAL(entirePipelineFinished()), this, SLOT(showPipelineFinishedLogMessage()));
    connect(scene, SIGNAL(entirePipelineFinished()), this, SLOT(updateMenu()));
    connect(scene, SIGNAL(pipelineExecutionFailed()), this, SLOT(updateMenu()));

    QRectF scene_rect = scene->itemsBoundingRect();
    tw->fitInView(scene_rect, Qt::KeepAspectRatio);
    tw->scale(0.75, 0.75);
    scene->setSceneRect(tw->mapToScene(tw->rect()).boundingRect());

    QRectF items_rect = scene->itemsBoundingRect();
    QRectF new_scene_rect = items_rect.united(tw->mapToScene(tw->rect()).boundingRect());
    qreal top_left_x = new_scene_rect.topLeft().x();
    qreal top_left_y = new_scene_rect.topLeft().y();
    qreal bottom_right_x = new_scene_rect.bottomRight().x();
    qreal bottom_right_y = new_scene_rect.bottomRight().y();
    qreal width = new_scene_rect.width();
    qreal height = new_scene_rect.height();
    new_scene_rect.setTopLeft(QPointF(top_left_x - width / 2.0, top_left_y - height / 2.0));
    new_scene_rect.setBottomRight(QPointF(bottom_right_x + width / 2.0, bottom_right_y + height / 2.0));
    scene->setSceneRect(new_scene_rect);

    desc_->blockSignals(true);
    desc_->setHtml(scene->getDescription());
    desc_->blockSignals(false);
  }

  void TOPPASBase::closeEvent(QCloseEvent* event)
  {
    QList<QMdiSubWindow*> all_windows = ws_->subWindowList();
    for (QMdiSubWindow* w : all_windows)
    {
      TOPPASWidget* widget = dynamic_cast<TOPPASWidget*>(w->widget());
      if (!widget) continue; // not a TOPPASWidget.. ignore it

      if (!widget->getScene()->saveIfChanged())
      { // user chose 'abort' in dialog
        event->ignore();
        return;
      }
    }
    event->accept();
    QSettings settings("OpenMS", "TOPPAS");
    settings.setValue("geometry", saveGeometry());
    settings.setValue("windowState", saveState());
  }

  void TOPPASBase::showURL()
  {
    QString target = dynamic_cast<QAction*>(sender())->data().toString();
    GUIHelpers::openURL(target);
  }

  TOPPASWidget* TOPPASBase::window_(int id) const
  {
    return dynamic_cast<TOPPASWidget*>(ws_->getWidget(id));
  }

  TOPPASWidget* TOPPASBase::activeSubWindow_() const
  {
    if (ws_ == nullptr || ws_->currentSubWindow() == nullptr)
    {
      return nullptr;
    }

    return dynamic_cast<TOPPASWidget*>(ws_->currentSubWindow()->widget());
  }

  void TOPPASBase::closeByTab(int id)
  {
    TOPPASWidget* window = window_(id);
    if (window)
    {
      // try to close the window.. the user might say no
      // The Dtor of 'window' will take care of removing it from the tabbar
      if (window->close()) updateMenu();
    }
  }

  void TOPPASBase::focusByTab(int id)
  {
    TOPPASWidget* window = window_(id);
    if (window)
    {
      //std::cerr << "tab changed...\n";
      desc_->blockSignals(true);
      desc_->setHtml(window->getScene()->getDescription());
      desc_->blockSignals(false);
      window->setFocus();
    }
    else
    {
      desc_->blockSignals(true);
      desc_->setHtml("");
      desc_->blockSignals(false);
    }
  }

  void TOPPASBase::closeFile()
  {
    if (ws_ != nullptr && ws_->currentSubWindow() != nullptr)
    {
      ws_->currentSubWindow()->close();
    }
    updateMenu();
  }

  void TOPPASBase::showStatusMessage(const string& msg, OpenMS::UInt time)
  {
    if (time == 0)
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

  void TOPPASBase::updateTabBar(QMdiSubWindow* w)
  {
    if (w)
    {
      TOPPASWidget* tw = dynamic_cast<TOPPASWidget*>(w->widget());
      if (tw)
      {
        Int window_id = tw->getWindowId();
        tab_bar_->show(window_id);
      }
    }
  }

  void TOPPASBase::loadPreferences(String filename)
  {
    //compose default ini file path
    String default_ini_file = String(QDir::homePath()) + "/.TOPPAS.ini";

    if (filename.empty())
    {
      filename = default_ini_file;
    }

    //load preferences, if file exists
    if (File::exists(filename))
    {
      bool error = false;
      Param tmp;
      ParamXMLFile paramFile;
      try // the file might be corrupt
      {
        paramFile.load(filename, tmp);
      }
      catch (...)
      {
        error = true;
      }

      //apply preferences if they are of the current TOPPAS version
      if (!error && tmp.exists("preferences:version") && tmp.getValue("preferences:version").toString() == VersionInfo::getVersion())
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
    param_.setValue("PreferencesFile", filename);

    // set the recent files
    recent_files_menu_.setFromParam(param_.copy("preferences:RecentFiles"));
  }

  void TOPPASBase::savePreferences()
  {
    // replace recent files
    param_.removeAll("preferences:RecentFiles");
    param_.insert("preferences:RecentFiles:", recent_files_menu_.getAsParam());

    //set version
    param_.setValue("preferences:version", VersionInfo::getVersion());

    Param save_param = param_.copy("preferences:");

    try
    {
      ParamXMLFile paramFile;
      // TODO: if closing multiple TOPPAS instances simultaneously, we might write to this file concurrently
      //       thus destroying its integrity. Think about using boost filelocks
      //       and also implement in TOPPView (and other GUI's which write to user directory)
      paramFile.store(string(param_.getValue("PreferencesFile")), save_param);
    }
    catch (Exception::UnableToCreateFile& /*e*/)
    {
      cerr << "Unable to create INI File: '" << string(param_.getValue("PreferencesFile")) << "'" << endl;
    }
  }

  void TOPPASBase::showAboutDialog()
  {
    QApplicationTOPP::showAboutDialog(this, "TOPPAS");
  }

  void TOPPASBase::updateMenu()
  {
    TOPPASWidget* tw = activeSubWindow_();
    TOPPASScene* ts = nullptr;
    if (tw)
    {
      ts = tw->getScene();
    }

    QList<QAction*> actions = this->findChildren<QAction*>("");
    for (int i = 0; i < actions.count(); ++i)
    {
      QString text = actions[i]->text();

      if (text == "&Run (F5)")
      {
        bool show = false;
        if (ts && !(ts->isPipelineRunning()))
        {
          show = true;
        }
        actions[i]->setEnabled(show);
      }
      else if (text == "&Abort")
      {
        bool show = false;
        if (ts && ts->isPipelineRunning())
        {
          show = true;
        }
        actions[i]->setEnabled(show);
      }
      else if (text == "&Include")
      {
        bool show = ts;
        actions[i]->setEnabled(show);
      }
      else if (text == "&Load resource file")
      {
        bool show = ts;
        actions[i]->setEnabled(show);
      }
      else if (text == "Save &resource file")
      {
        bool show = ts;
        actions[i]->setEnabled(show);
      }
      else if (text == "&Save")
      {
        bool show = ts && ts->wasChanged();
        actions[i]->setEnabled(show);
      }
      else if (text == "Refresh &parameters")
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
        title = asterisk_shown ? title.right(title.size() - 1) : QString("*") + title;
        tw->setWindowTitle(title);
        tab_bar_->setTabText(title);
      }
    }
  }

  void TOPPASBase::keyPressEvent(QKeyEvent* e)
  {
    if (e->key() == Qt::Key_F5)
    {
      TOPPASWidget* tw = activeSubWindow_();
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
    if (param_.getValue("preferences:default_path_current") != "true")
      return;

    //reset
    current_path_ = param_.getValue("preferences:default_path").toString();

    //update if the current layer has a path associated TODO
    //if (activeCanvas_() && activeCanvas_()->getLayerCount()!=0 && activeCanvas_()->getCurrentLayer().filename!="")
    //{
    // current_path_ = File::path(activeCanvas_()->getCurrentLayer().filename);
    //}
  }

  void TOPPASBase::insertNewVertex_(double x, double y, QTreeWidgetItem* item)
  {
    if (!activeSubWindow_() || !activeSubWindow_()->getScene() || !tools_tree_view_)
    {
      return;
    }

    TOPPASScene* scene = activeSubWindow_()->getScene();
    QTreeWidgetItem* current_tool = item ? item : tools_tree_view_->currentItem();
    String tool_name = String(current_tool->text(0));
    TOPPASVertex* tv = nullptr;

    if (tool_name == "<Input files>")
    {
      tv = new TOPPASInputFileListVertex();
    }
    else if (tool_name == "<Output files>")
    {
      tv = new TOPPASOutputFileListVertex();
      TOPPASOutputFileListVertex* oflv = dynamic_cast<TOPPASOutputFileListVertex*>(tv);
      connect(tv, SIGNAL(outputFileWritten(const String &)), this, SLOT(outputVertexFinished(const String &)));
      scene->connectOutputVertexSignals((TOPPASOutputVertex*)oflv);
    }
    else if (tool_name == "<Output folder>")
    {
      tv = new TOPPASOutputFolderVertex();
      TOPPASOutputFolderVertex* oflv = dynamic_cast<TOPPASOutputFolderVertex*>(tv);
      connect(tv, SIGNAL(outputFileWritten(const String&)), this, SLOT(outputVertexFinished(const String&)));
      scene->connectOutputVertexSignals((TOPPASOutputVertex*)oflv);
    }
    else if (tool_name == "<Merger>")
    {
      tv = new TOPPASMergerVertex(true);
      connect(tv, SIGNAL(mergeFailed(const QString)), this, SLOT(updateTOPPOutputLog(const QString &)));
    }
    else if (tool_name == "<Collector>")
    {
      tv = new TOPPASMergerVertex(false);
      connect(tv, SIGNAL(mergeFailed(const QString)), this, SLOT(updateTOPPOutputLog(const QString &)));
    }
    else if (tool_name == "<Splitter>")
    {
      tv = new TOPPASSplitterVertex();
    }
    else // node is a TOPP tool
    {
      if (current_tool->childCount() > 0)
      {
        // category or tool name with types is selected (instead of a concrete type)
        return;
      }
      String tool_type;
      if (current_tool->parent() != nullptr && current_tool->parent()->parent() != nullptr)
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
      TOPPASToolVertex* ttv = dynamic_cast<TOPPASToolVertex*>(tv);

      // check if tool init was successful (i.e. tool was found); TODO: only populate Tool list with available tools so we do not need to check?!
      if (!ttv->isToolReady())
      {
        delete ttv;
        return;
      }

      connect(ttv, SIGNAL(toolStarted()), this, SLOT(toolStarted()));
      connect(ttv, SIGNAL(toolFinished()), this, SLOT(toolFinished()));
      connect(ttv, SIGNAL(toolCrashed()), this, SLOT(toolCrashed()));
      connect(ttv, SIGNAL(toolFailed()), this, SLOT(toolFailed()));
      // already done in ToppasScene:
      //connect (ttv, SIGNAL(toppOutputReady(const QString&)), this, SLOT(updateTOPPOutputLog(const QString&)));

      scene->connectToolVertexSignals(ttv);
    }

    scene->connectVertexSignals(tv);
    scene->addVertex(tv);
    tv->setPos(x, y);
    tv->setZValue(z_value_);
    z_value_ += 0.000001;
    scene->topoSort(false);
    scene->setChanged(true);
  }

  void TOPPASBase::runPipeline()
  {
    TOPPASWidget* w = activeSubWindow_();
    if (w)
    {
      w->getScene()->runPipeline();
    }
  }

  void TOPPASBase::abortPipeline()
  {
    TOPPASWidget* w = activeSubWindow_();
    if (w)
    {
      w->getScene()->abortPipeline();
    }
    updateMenu();
  }

  void TOPPASBase::toolStarted()
  {
    TOPPASToolVertex* tv = dynamic_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " of node #" + String(tv->getTopoNr()) + " started. Processing ...";

      log_->appendNewHeader(LogWindow::LogState::NOTICE, text, "");
    }
    updateMenu();
  }

  void TOPPASBase::toolFinished()
  {
    TOPPASToolVertex* tv = dynamic_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " finished!";

      log_->appendNewHeader(LogWindow::LogState::NOTICE, text, "");
    }
    updateMenu();
  }

  void TOPPASBase::toolCrashed()
  {
    TOPPASToolVertex* tv = dynamic_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " crashed!";

      log_->appendNewHeader(LogWindow::LogState::CRITICAL, text, "");
    }
    updateMenu();
  }

  void TOPPASBase::toolFailed()
  {
    TOPPASToolVertex* tv = dynamic_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " failed!";

      log_->appendNewHeader(LogWindow::LogState::CRITICAL, text, "");
    }
    updateMenu();
  }

  void TOPPASBase::outputVertexFinished(const String& file)
  {
    String text = "Output file '" + file + "' written.";
    log_->appendNewHeader(LogWindow::LogState::NOTICE, text, "");
  }

  void TOPPASBase::updateTOPPOutputLog(const QString& out)
  {
    const QString& text = out; // shortened version for now (if we reintroduce simultaneous tool execution,
                        // we need to rethink this (probably only trigger this slot when tool 100% finished)


    //show log if there is output
    dynamic_cast<QWidget*>(log_->parent())->show();

    //update log_
    log_->moveCursor(QTextCursor::End, QTextCursor::MoveAnchor); // move cursor to end, since text is inserted at cursor
    log_->insertPlainText(text);
  }

  void TOPPASBase::showPipelineFinishedLogMessage()
  {
    log_->appendNewHeader(LogWindow::LogState::NOTICE, "Entire pipeline execution finished!", "");
  }

  void TOPPASBase::insertNewVertexInCenter_(QTreeWidgetItem* item)
  {
    if (!activeSubWindow_() || !activeSubWindow_()->getScene() || !tools_tree_view_ || !tools_tree_view_->currentItem())
    {
      return;
    }

    QPointF insert_pos = activeSubWindow_()->mapToScene(QPoint((activeSubWindow_()->width() / 2.0) + (qreal)(5 * node_offset_), (activeSubWindow_()->height() / 2.0) + (qreal)(5 * node_offset_)));
    insertNewVertex_(insert_pos.x(), insert_pos.y(), item);
    node_offset_ = (node_offset_ + 1) % 10;
  }

  void TOPPASBase::saveToClipboard(TOPPASScene* scene)
  {
    if (clipboard_scene_ != nullptr)
    {
      delete clipboard_scene_;
      clipboard_scene_ = nullptr;
    }
    clipboard_scene_ = scene;
  }

  void TOPPASBase::sendClipboardContent()
  {
    TOPPASScene* sndr = dynamic_cast<TOPPASScene*>(QObject::sender());
    if (sndr != nullptr)
    {
      sndr->setClipboard(clipboard_scene_);
    }
  }

  void TOPPASBase::refreshParameters()
  {
    TOPPASWidget* w = activeSubWindow_();
    QString file_name = TOPPASBase::refreshPipelineParameters(w, current_path_.toQString());
    if (file_name != "")
    {
      tab_bar_->setTabText(File::basename(file_name).toQString());
    }
  }

  // static
  QString TOPPASBase::refreshPipelineParameters(TOPPASWidget* tw, QString current_path)
  {
    TOPPASScene* ts = nullptr;
    if (tw)
    {
      ts = tw->getScene();
    }
    if (!ts)
    {
      return "";
    }

    TOPPASScene::RefreshStatus st = ts->refreshParameters();
    if (st == TOPPASScene::ST_REFRESH_NOCHANGE)
    {
      QMessageBox::information(tw, tr("Nothing to be done"),
                               tr("The parameters of the tools used in this workflow have not changed."));
      return "";
    }

    ts->setChanged(true);
    ts->updateEdgeColors();
    if (st == TOPPASScene::ST_REFRESH_CHANGEINVALID)
    {
      QMessageBox::information(tw, "Parameters updated!",
                               "The resulting pipeline is now invalid. Probably some input or output parameters were removed or added. Please repair!",
                               QMessageBox::Ok);
      return "";
    }
    else if (st == TOPPASScene::ST_REFRESH_REMAINSINVALID)
    {
      QMessageBox::information(tw, "Parameters updated!",
                               "The resulting pipeline remains invalid (not runnable). Maybe some input files or even edges are missing. Please repair!",
                               QMessageBox::Ok);
      return "";
    }

    int ret = QMessageBox::information(tw, "Parameters updated!",
                                       "The parameters of some tools in this workflow have changed. Do you want to save these changes now?",
                                       QMessageBox::Save | QMessageBox::Cancel);
    if (ret == QMessageBox::Save)
    {
      QString file_name = TOPPASBase::savePipelineAs(tw, std::move(current_path));
      return file_name;
    }

    return "";
  }

  void TOPPASBase::openFilesInTOPPView(QStringList files)
  {
    if (files.empty()) return;
    
    if (files.size() > 1)
    {
      // ask user how to open multiple files
      QMessageBox msgBox(
        QMessageBox::Question,
        tr("Open files with overlay?"),
        tr("How do you want to open the output files?"),
        QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
      msgBox.setButtonText(QMessageBox::Yes, tr("&Single Tab - Overlay"));
      msgBox.setButtonText(QMessageBox::No, tr("&Separate tabs"));
      int ret = msgBox.exec();
      if (ret == QMessageBox::Cancel) return; // Escape was pressed
      if (ret == QMessageBox::Yes)
      { // put a '+' in between the files (TOPPView's command line will interpret this as overlay)
        files = files.join("#SpLiT_sTrInG#+#SpLiT_sTrInG#").split("#SpLiT_sTrInG#", Qt::SkipEmptyParts);
      }
    }
    
    if (!GUIHelpers::startTOPPView(files))
    {
      QMessageBox::warning(this, "Could not start TOPPView", "TOPPView failed to start. Please see the commandline for details.");
    }

  }

  void TOPPASBase::openToppasFile(const QString& filename)
  {
    addTOPPASFile(String(filename));
  }

} //namespace OpenMS

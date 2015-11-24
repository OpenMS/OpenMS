// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <cstdio>
#include <cstdlib>

#include <OpenMS/VISUAL/APPLICATIONS/TOPPASBase.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASLogWindow.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASResources.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASSplitterVertex.h>
#include <OpenMS/VISUAL/TOPPASTabBar.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASWidget.h>

//Qt
#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QMap>
#include <QtCore/QSet>
#include <QtCore/QUrl>
#include <QtGui/QApplication>
#include <QtGui/QCheckBox>
#include <QtGui/QCloseEvent>
#include <QtGui/QDesktopServices>
#include <QtGui/QDesktopWidget>
#include <QtGui/QDockWidget>
#include <QtGui/QFileDialog>
#include <QtGui/QInputDialog>
#include <QtGui/QLabel>
#include <QtGui/QListWidget>
#include <QtGui/QListWidgetItem>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QMessageBox>
#include <QtGui/QSplashScreen>
#include <QtGui/QStatusBar>
#include <QtGui/QTextEdit>
#include <QtGui/QToolBar>
#include <QtGui/QToolButton>
#include <QtGui/QTreeWidget>
#include <QtGui/QTreeWidgetItem>
#include <QtGui/QToolTip>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWhatsThis>

#include <QWebView>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QTextStream>
#include <QSvgGenerator>
#include <QNetworkProxyFactory>
#include <QNetworkProxy>
#include <QTextCodec>

using namespace std;

namespace OpenMS
{
  using namespace Internal;

  int TOPPASBase::node_offset_ = 0;
  qreal TOPPASBase::z_value_ = 42.0;

  TOPPASBase::TOPPASBase(QWidget* parent) :
    QMainWindow(parent),
    DefaultParamHandler("TOPPASBase"),
    clipboard_scene_(0)
  {
#if  defined(__APPLE__)
    // we do not want to load plugins as this leads to serious problems
    // when shipping on mac os x
    QApplication::setLibraryPaths(QStringList());
#endif

    setWindowTitle("TOPPAS");
    setWindowIcon(QIcon(":/TOPPAS.png"));

    // ensure correct encoding of paths
    QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF-8"));

    //prevents errors caused by too small width,height values
    setMinimumSize(400, 400);

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
    tab_bar_->addTab("dummy", 1336);
    tab_bar_->setMinimumSize(tab_bar_->sizeHint());
    tab_bar_->removeId(1336);
    //connect slots and signals for selecting spectra
    connect(tab_bar_, SIGNAL(currentIdChanged(int)), this, SLOT(focusByTab(int)));
    connect(tab_bar_, SIGNAL(aboutToCloseId(int)), this, SLOT(closeByTab(int)));

    box_layout->addWidget(tab_bar_);
    ws_ = new QWorkspace(dummy);
    connect(ws_, SIGNAL(windowActivated(QWidget*)), this, SLOT(updateTabBar(QWidget*)));
    connect(ws_, SIGNAL(windowActivated(QWidget*)), this, SLOT(updateMenu()));

    box_layout->addWidget(ws_);

    //################## MENUS #################
    // File menu
    QMenu* file = new QMenu("&File", this);
    menuBar()->addMenu(file);
    file->addAction("&New", this, SLOT(newPipeline()), Qt::CTRL + Qt::Key_N);
    file->addAction("&Open", this, SLOT(openFileDialog()), Qt::CTRL + Qt::Key_O);
    file->addAction("Open &example file", this, SLOT(openExampleDialog()), Qt::CTRL + Qt::Key_E);
    file->addAction("&Include", this, SLOT(includePipeline()), Qt::CTRL + Qt::Key_I);
    file->addAction("Online &Repository", this, SLOT(openOnlinePipelineRepository()), Qt::CTRL + Qt::Key_R);
    file->addAction("&Save", this, SLOT(savePipeline()), Qt::CTRL + Qt::Key_S);
    file->addAction("Save &As", this, SLOT(saveCurrentPipelineAs()), Qt::CTRL + Qt::SHIFT + Qt::Key_S);
    file->addAction("E&xport as image", this, SLOT(exportAsImage()));
    file->addAction("Refresh &parameters", this, SLOT(refreshParameters()), Qt::CTRL + Qt::SHIFT + Qt::Key_P);
    file->addAction("&Close", this, SLOT(closeFile()), Qt::CTRL + Qt::Key_W);
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
    defaults_.setValidStrings("preferences:default_path_current", ListUtils::create<String>("true,false"));
    defaults_.setValue("preferences:version", "none", "OpenMS version, used to check if the TOPPAS.ini is up-to-date");

    defaultsToParam_();

    //load param file
    loadPreferences();

    //################## Dock widgets #################
    //TOPP tools window
    QDockWidget* topp_tools_bar = new QDockWidget("TOPP", this);
    addDockWidget(Qt::LeftDockWidgetArea, topp_tools_bar);
    tools_tree_view_ = createTOPPToolsTreeWidget(topp_tools_bar);
    topp_tools_bar->setWidget(tools_tree_view_);
    connect(tools_tree_view_, SIGNAL(itemDoubleClicked(QTreeWidgetItem*, int)), this, SLOT(insertNewVertexInCenter_(QTreeWidgetItem*)));
    windows->addAction(topp_tools_bar->toggleViewAction());

    //log window
    QDockWidget* log_bar = new QDockWidget("Log", this);
    addDockWidget(Qt::BottomDockWidgetArea, log_bar);
    log_ = new TOPPASLogWindow(log_bar);
    log_->setReadOnly(true);
    log_->setMaxLength(1e7); // limit to 10 mio characters, and trim to 5 mio upon reaching this limit
    log_bar->setWidget(log_);
    log_bar->hide();
    //windows->addAction("&Show log window",log_bar,SLOT(show()));
    windows->addAction(log_bar->toggleViewAction());

    //workflow description window
    QDockWidget* description_bar = new QDockWidget("Workflow Description", this);
    addDockWidget(Qt::RightDockWidgetArea, description_bar);
    desc_ = new QTextEdit(description_bar);
    desc_->setTextColor(Qt::black);
    desc_->setText("... put your workflow description here ...");
    desc_->setTextColor(Qt::black);
    desc_->document()->setDefaultFont(QFont("Arial", 12));
    description_bar->setWidget(desc_);
    //windows->addAction("&Show log window",log_bar,SLOT(show()));
    windows->addAction(description_bar->toggleViewAction());
    connect(desc_, SIGNAL(textChanged()), this, SLOT(descriptionUpdated_()));

    // set current path
    current_path_ = param_.getValue("preferences:default_path");

    // set & create temporary path -- make sure its a new subdirectory, as it will be deleted later
    QString new_tmp_dir = File::getUniqueName().toQString();
    QDir qd(File::getTempDirectory().toQString());
    qd.mkdir(new_tmp_dir);
    qd.cd(new_tmp_dir);
    tmp_path_ = qd.absolutePath();

    // online browser
    webview_ = new QWebView(parent);
    webview_->page()->setLinkDelegationPolicy(QWebPage::DelegateAllLinks); // now linkClicked() is emitted

    connect((webview_->page()), SIGNAL(linkClicked(const QUrl &)), this, SLOT(downloadTOPPASfromHomepage_(const QUrl &)));

    network_manager_ = new QNetworkAccessManager(this);
    connect(network_manager_, SIGNAL(finished(QNetworkReply*)), this, SLOT(toppasFileDownloaded_(QNetworkReply*)));

    // update the menu
    updateMenu(); 
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
    if (!activeWindow_() || !activeWindow_()->getScene())
    {
      return;
    }
    //std::cerr << "changed to '" << String(desc_->toHtml()) << "'\n";
    activeWindow_()->getScene()->setChanged(true);
    activeWindow_()->getScene()->setDescription(desc_->toHtml());
  }

  void TOPPASBase::toppasFileDownloaded_(QNetworkReply* r)
  {
    r->deleteLater();
    if (r->error() != QNetworkReply::NoError)
    {
      showLogMessage_(LS_ERROR, "Download failed", "Error '" + r->errorString() + "' while downloading TOPPAS file: '" + r->url().toString() + "'");
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
      LOG_WARN << "The URL format of downloads from the TOPPAS Online-Repository has changed. Please notify developers!";
    }
    QString filename = QFileDialog::getSaveFileName(this, "Where to save the TOPPAS file?", this->current_path_.toQString() + "/" + proposed_filename, tr("TOPPAS (*.toppas)"));

    // check if the user clicked cancel, to avoid saving .toppas somewhere
    if (String(filename).trim().empty())
    {
      showLogMessage_(LS_NOTICE, "Download succeeded, but saving aborted by user!", "");
      return;
    }

    if (!filename.endsWith(".toppas", Qt::CaseInsensitive))
    {
      filename += ".toppas";
    }

    QFile file(filename);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
      showLogMessage_(LS_NOTICE, "Download succeeded. Cannot save the file. Try again with another filename and/or location!", "");
      return;
    }

    QTextStream out(&file);
    out << data;
    file.close();

    this->addTOPPASFile(filename);
    showLogMessage_(LS_NOTICE, "File successfully saved to '" + filename + "'.", "");
  }
  
  void TOPPASBase::TOPPASreadyRead()
  {
    QNetworkReply::NetworkError ne = network_reply_->error();
    qint64 ba = network_reply_->bytesAvailable();
    LOG_DEBUG << "Error code (QNetworkReply::NetworkError): " << ne << "  bytes available: " << ba << std::endl;
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

      showLogMessage_(LS_NOTICE, "Downloading file '" + url.toString() + "'. You will be notified once the download finished.", "");
      webview_->close();
    }
    else
    {
      QMessageBox::warning(this, tr("Error"), tr("You can only click '.toppas' files on this page. No navigation is allowed!\n"));
      webview_->setFocus();
      webview_->activateWindow();
    }
  }

  void TOPPASBase::openOnlinePipelineRepository()
  {
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
    item = new QTreeWidgetItem((QTreeWidget*)0);
    item->setText(0, "<Collector>");
    tools_tree_view->addTopLevelItem(item);
    item = new QTreeWidgetItem((QTreeWidget*)0);
    item->setText(0, "<Splitter>");
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
      if (it->second.category.trim() == "")
        it->second.category = "Unassigned";
    }

    QSet<QString> category_set;
    for (ToolListType::ConstIterator it = tools_list.begin(); it != tools_list.end(); ++it)
    {
      category_set << String(it->second.category).toQString();
    }
    QStringList category_list = category_set.toList();
    qSort(category_list);
    Map<QString, QTreeWidgetItem*> category_map;

    foreach(const QString &category, category_list)
    {
      item = new QTreeWidgetItem((QTreeWidget*)0);
      item->setText(0, category);
      tools_tree_view->addTopLevelItem(item);
      category_map[category] = item;
    }

    for (ToolListType::iterator it = tools_list.begin(); it != tools_list.end(); ++it)
    {
      item = new QTreeWidgetItem(category_map[it->second.category.toQString()]);
      item->setText(0, it->first.toQString());
      QTreeWidgetItem* parent_item = item;
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
    if (file_name == "") return;

    if (!file_name.toQString().endsWith(".toppas", Qt::CaseInsensitive))
    {
      LOG_ERROR << "The file '" << file_name << "' is not a .toppas file" << std::endl;
      return;
    }

    TOPPASScene* scene = 0;
    if (in_new_window)
    {
      if (activeWindow_())
      {
        TOPPASWidget* uninitialized_window = window_(IDINITIALUNTITLED);
        if (uninitialized_window && !uninitialized_window->getScene()->wasChanged())
          closeByTab(IDINITIALUNTITLED);
      }
      TOPPASWidget* tw = new TOPPASWidget(Param(), ws_, tmp_path_);
      scene = tw->getScene();
      scene->load(file_name); // first load WF, including description etc
      showAsWindow_(tw, File::basename(file_name)); // show it
    }
    else
    {
      if (!activeWindow_()) return;

      TOPPASScene* tmp_scene = new TOPPASScene(0, this->tmp_path_.toQString(), false);
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
        connect(tv, SIGNAL(toolStarted()), this, SLOT(toolStarted()));
        connect(tv, SIGNAL(toolFinished()), this, SLOT(toolFinished()));
        connect(tv, SIGNAL(toolCrashed()), this, SLOT(toolCrashed()));
        connect(tv, SIGNAL(toolFailed()), this, SLOT(toolFailed()));
        connect(tv, SIGNAL(toolFailed(const QString &)), this, SLOT(updateTOPPOutputLog(const QString &)));
        // already done in ToppasScene:
        //connect (tv, SIGNAL(toppOutputReady(const QString&)), this, SLOT(updateTOPPOutputLog(const QString&)));
        continue;
      }

      TOPPASMergerVertex* tmv = qobject_cast<TOPPASMergerVertex*>(*it);
      if (tmv)
      {
        connect(tmv, SIGNAL(mergeFailed(const QString)), this, SLOT(updateTOPPOutputLog(const QString &)));
        continue;
      }

      TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(*it);
      if (oflv)
      {
        connect(oflv, SIGNAL(outputFileWritten(const String &)), this, SLOT(outputVertexFinished(const String &)));
        continue;
      }
    }
  }

  void TOPPASBase::newPipeline(const int id)
  {
    TOPPASWidget* tw = new TOPPASWidget(Param(), ws_, tmp_path_);
    showAsWindow_(tw, "(Untitled)", id);
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
      if (!w->getScene()->store(file_name))
      {
        QMessageBox::warning(NULL, tr("Error"),
                             tr("Unable to save current pipeline. Possible reason: Invalid edges due to parameter refresh."));
      }
      QString caption = File::basename(file_name).toQString();
      w->setWindowTitle(caption);
    }
    return file_name;
  }

  void TOPPASBase::exportAsImage()
  {
    TOPPASWidget* w = activeWindow_();
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

  void TOPPASBase::showAsWindow_(TOPPASWidget* tw, const String& caption, const int special_id)
  {
    ws_->addWindow(tw);
    connect(tw, SIGNAL(sendStatusMessage(std::string, OpenMS::UInt)), this, SLOT(showStatusMessage(std::string, OpenMS::UInt)));
    connect(tw, SIGNAL(sendCursorStatus(double, double)), this, SLOT(showCursorStatus(double, double)));
    connect(tw, SIGNAL(toolDroppedOnWidget(double, double)), this, SLOT(insertNewVertex_(double, double)));
    connect(tw, SIGNAL(pipelineDroppedOnWidget(const String &, bool)), this, SLOT(addTOPPASFile(const String &, bool)));
    tw->setWindowTitle(caption.toQString());

    //add tab with id
    static int window_counter = 1337;
    ++window_counter;

    // use special_id if given (for first untitled tab), otherwise the running window_counter
    int local_counter = special_id == -1 ? window_counter : special_id;
    tw->setWindowId(local_counter);

    tab_bar_->addTab(caption.toQString(), tw->getWindowId());

    //connect slots and signals for removing the widget from the bar, when it is closed
    //- through the menu entry
    //- through the tab bar
    //- through the MDI close button
    connect(tw, SIGNAL(aboutToBeDestroyed(int)), tab_bar_, SLOT(removeId(int)));

    tab_bar_->setCurrentId(tw->getWindowId());

    //show first window maximized (only visible windows are in the list)
    if (ws_->windowList().count() == 0)
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
    bool close = true;
    QList<QWidget*> all_windows = ws_->windowList();
    foreach(QWidget * w, all_windows)
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
    QString target = qobject_cast<QAction*>(sender())->data().toString();
    GUIHelpers::openURL(target);
  }

  TOPPASWidget* TOPPASBase::window_(int id) const
  {
    //cout << "Looking for tab with id: " << id << endl;
    QList<QWidget*> windows = ws_->windowList();
    for (int i = 0; i < windows.size(); ++i)
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
    if (!ws_->activeWindow())
      return 0;

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
    ws_->activeWindow()->close();
    updateMenu();
  }

  void TOPPASBase::showStatusMessage(string msg, OpenMS::UInt time)
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

    if (filename == "")
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
  }

  void TOPPASBase::savePreferences()
  {
    //set version
    param_.setValue("preferences:version", VersionInfo::getVersion());

    Param save_param = param_.copy("preferences:");

    try
    {
      ParamXMLFile paramFile;
      // TODO: if closing multiple TOPPAS instances simultaneously, we might write to this file concurrently
      //       thus destroying its integrity. Think about using boost filelocks
      //       see OpenMS/METADATA/DocumentIDTagger.h for example
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
    //dialog and grid layout
    QDialog* dlg = new QDialog(this);
    QGridLayout* grid = new QGridLayout(dlg);
    dlg->setWindowTitle("About TOPPAS");

    QLabel* label = new QLabel(dlg);
    label->setPixmap(QPixmap(":/TOPP_about.png"));
    grid->addWidget(label, 0, 0);

    //text
    QString text = QString("<BR>"
                           "<FONT size=+3>TOPPAS</font><BR>"
                           "<BR>"
                           "Version: %1%2<BR>"
                           "<BR>"
                           "OpenMS and TOPP is free software available under the<BR>"
                           "BSD 3-Clause Licence (BSD-new)<BR>"
                           "<BR>"
                           "<BR>"
                           "<BR>"
                           "<BR>"
                           "<BR>"
                           "Any published work based on TOPP and OpenMS shall cite these papers:<BR>"
                           "Sturm et al., BMC Bioinformatics (2008), 9, 163<BR>"
                           "Kohlbacher et al., Bioinformatics (2007), 23:e191-e197<BR>"
                           ).arg(VersionInfo::getVersion().toQString()
                           ).arg( // if we have a revision, embed it also into the shown version number
                              VersionInfo::getRevision() != "" ? QString(" (") + VersionInfo::getRevision().toQString() + ")" : "");
    
    QLabel* text_label = new QLabel(text, dlg);
    grid->addWidget(text_label, 0, 1, Qt::AlignTop | Qt::AlignLeft);

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
        tab_bar_->setTabText(tab_bar_->currentIndex(), title);
      }
    }
  }

  void TOPPASBase::showLogMessage_(TOPPASBase::LogState state, const String& heading, const String& body)
  {
    //Compose current time string
    DateTime d = DateTime::now();

    String state_string;
    switch (state)
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
    qobject_cast<QWidget*>(log_->parent())->show();
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
    if (param_.getValue("preferences:default_path_current") != "true")
      return;

    //reset
    current_path_ = param_.getValue("preferences:default_path");

    //update if the current layer has a path associated TODO
    //if (activeCanvas_() && activeCanvas_()->getLayerCount()!=0 && activeCanvas_()->getCurrentLayer().filename!="")
    //{
    // current_path_ = File::path(activeCanvas_()->getCurrentLayer().filename);
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
      connect(oflv, SIGNAL(outputFileWritten(const String &)), this, SLOT(outputVertexFinished(const String &)));
      scene->connectOutputVertexSignals(oflv);
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

      tv = new TOPPASToolVertex(tool_name, tool_type);
      TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(tv);

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
        text += " (" + type + ")";
      }
      text += " of node #" + String(tv->getTopoNr()) + " started. Processing ...";

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
        text += " (" + type + ")";
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
        text += " (" + type + ")";
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
        text += " (" + type + ")";
      }
      text += " failed!";

      showLogMessage_(LS_ERROR, text, "");
    }
    updateMenu();
  }

  void TOPPASBase::outputVertexFinished(const String& file)
  {
    String text = "Output file '" + file + "' written.";
    showLogMessage_(LS_NOTICE, text, "");
  }

  void TOPPASBase::updateTOPPOutputLog(const QString& out)
  {
    QString text = out; // shortened version for now (if we reintroduce simultaneous tool execution,
                        // we need to rethink this (probably only trigger this slot when tool 100% finished)


    //show log if there is output
    qobject_cast<QWidget*>(log_->parent())->show();

    //update log_
    log_->moveCursor(QTextCursor::End, QTextCursor::MoveAnchor); // move cursor to end, since text is inserted at cursor
    log_->insertPlainText(text);
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

    QPointF insert_pos = activeWindow_()->mapToScene(QPoint((activeWindow_()->width() / 2.0) + (qreal)(5 * node_offset_), (activeWindow_()->height() / 2.0) + (qreal)(5 * node_offset_)));
    insertNewVertex_(insert_pos.x(), insert_pos.y(), item);
    node_offset_ = (node_offset_ + 1) % 10;
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
      QString file_name = TOPPASBase::savePipelineAs(tw, current_path);
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
      {
        files = files.join("#SpLiT_sTrInG#+#SpLiT_sTrInG#").split("#SpLiT_sTrInG#", QString::SkipEmptyParts);
      }
    }
    
    GUIHelpers::startTOPPView(files);

  }

  void TOPPASBase::openToppasFile(QString filename)
  {
    addTOPPASFile(String(filename));
  }

} //namespace OpenMS

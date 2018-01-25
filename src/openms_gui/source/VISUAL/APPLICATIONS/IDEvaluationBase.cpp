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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <cstdio>
#include <cstdlib>

#include <OpenMS/VISUAL/APPLICATIONS/IDEvaluationBase.h>

#include <OpenMS/CONCEPT/VersionInfo.h>

#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASWidget.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASTabBar.h>
#include <OpenMS/VISUAL/TOPPASResources.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>


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

#include <QWebView>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QTextStream>
#include <QSvgGenerator>
#include <QNetworkProxyFactory>
#include <QNetworkProxy>

using namespace std;

namespace OpenMS
{
  using namespace Internal;

  IDEvaluationBase::IDEvaluationBase(QWidget* parent) :
    QMainWindow(parent),
    DefaultParamHandler("IDEvaluationBase"),
    spec_1d_(nullptr)
  {
    for (double d = 0.0; d <= 1.0; d += (1.0) / 100)
    {
      q_value_thresholds_.push_back(d);
    }

    setWindowTitle("IDEvaluationBase");
    setWindowIcon(QIcon(":/TOPPAS.png"));

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

    ws_ = new QWorkspace(dummy);
    //connect(ws_,SIGNAL(windowActivated(QWidget*)),this,SLOT(updateMenu()));

    box_layout->addWidget(ws_);

    //################## MENUS #################
    // File menu
    QMenu* file = new QMenu("&File", this);
    menuBar()->addMenu(file);
    file->addAction("Add search result", this, SLOT(openFileDialog()), Qt::CTRL + Qt::Key_O);
    file->addAction("Save Image &As", this, SLOT(saveImageAs()), Qt::CTRL + Qt::Key_S);
    file->addSeparator();
    file->addAction("&Quit", qApp, SLOT(quit()));

    //Advanced menu
    //QMenu* advanced = new QMenu("&Advanced",this);
    //menuBar()->addMenu(advanced);
    //advanced->addAction("&Refresh definitions",this,SLOT(refreshDefinitions()), Qt::CTRL+Qt::Key_R);

    //Help menu
    QMenu* help = new QMenu("&Help", this);
    menuBar()->addMenu(help);
    QAction* action = help->addAction("OpenMS website", this, SLOT(showURL()));
    action->setData("http://www.OpenMS.de");
    help->addAction("&About", this, SLOT(showAboutDialog()));

    //create status bar
    message_label_ = new QLabel(statusBar());
    statusBar()->addWidget(message_label_, 1);


    //################## Dock widgets #################
    //TOPP tools window
    /*QDockWidget* plot_bar = new QDockWidget("Plot", this);
    addDockWidget(Qt::AllDockWidgetAreas, plot_bar);
    plot_bar->showMaximized();*/
    spec_1d_ =  new Spectrum1DWidget(Param(), this);
    spec_1d_->xAxis()->setLegend("q-value");
    // y axis is done in setIntensityMode();
    Param legend_on;
    legend_on.setValue("show_legend", "true", "Annotate each layer with its name on the canvas.");
    spec_1d_->canvas()->setParameters(legend_on);

    this->setCentralWidget(spec_1d_);

    //log window
    QDockWidget* log_bar = new QDockWidget("Log", this);
    addDockWidget(Qt::BottomDockWidgetArea, log_bar);
    log_ = new QTextEdit(log_bar);
    log_->setReadOnly(true);
    log_bar->setWidget(log_);
    log_bar->hide();

    //################## TOOLBARS #################
    //create toolbars and connect signals
    //--Basic tool bar for all views--
    tool_bar_ = addToolBar("Basic tool bar");

    //intensity modes
    intensity_button_group_ = new QButtonGroup(tool_bar_);
    intensity_button_group_->setExclusive(true);

    QToolButton* b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/lin.png"));
    b->setToolTip("PSM-Count: Normal");
    b->setShortcut(Qt::Key_N);
    b->setCheckable(true);
    b->setWhatsThis("PSM-Count: Normal<BR><BR>PSM-Count is displayed unmodified.<BR>(Hotkey: N)");
    intensity_button_group_->addButton(b, SpectrumCanvas::IM_NONE);
    tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/percentage.png"));
    b->setToolTip("PSM-Count: Percentage");
    b->setShortcut(Qt::Key_P);
    b->setCheckable(true);
    b->setWhatsThis("PSM-Count: Percentage<BR><BR>PSM-Count is displayed as a percentage of the layer"
                    " maximum PSM-Count. If only one layer is displayed this mode behaves like the"
                    " normal mode. If more than one layer is displayed PSM-Count are aligned."
                    "<BR>(Hotkey: P)");
    intensity_button_group_->addButton(b, SpectrumCanvas::IM_PERCENTAGE);
    tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/snap.png"));
    b->setToolTip("PSM-Count: Snap to maximum displayed PSM-Count");
    b->setShortcut(Qt::Key_S);
    b->setCheckable(true);
    b->setWhatsThis("PSM-Count: Snap to maximum displayed PSM-Count"
                    "<BR>(Hotkey: S)");
    intensity_button_group_->addButton(b, SpectrumCanvas::IM_SNAP);
    tool_bar_->addWidget(b);

    b = new QToolButton(tool_bar_);
    b->setIcon(QIcon(":/log.png"));
    b->setToolTip("PSM-Count: Use log scaling");
    b->setCheckable(true);
    b->setWhatsThis("PSM-Count: Logarithmic scaling of intensities for color calculation");
    intensity_button_group_->addButton(b, SpectrumCanvas::IM_LOG);
    tool_bar_->addWidget(b);

    connect(intensity_button_group_, SIGNAL(buttonClicked(int)), this, SLOT(setIntensityMode(int)));
    tool_bar_->addSeparator();

    //common buttons
    QAction* reset_zoom_button = tool_bar_->addAction(QIcon(":/reset_zoom.png"), "Reset Zoom", this, SLOT(resetZoom()));
    reset_zoom_button->setWhatsThis("Reset zoom: Zooms out as far as possible and resets the zoom history.<BR>(Hotkey: Backspace)");

    tool_bar_->show();


    //// parameters
    defaults_.setValue("image:height", 800, "Height of raster images (e.g., PNG).");
    defaults_.setValue("image:width", 1024, "Width of raster images (e.g., PNG).");
    FalseDiscoveryRate fdr;
    defaults_.insert("fdr:", fdr.getParameters());
    defaultsToParam_();
    //update the menu
    //updateMenu();
  }

  IDEvaluationBase::~IDEvaluationBase()
  {
  }

  QSize IDEvaluationBase::sizeHint() const
  {
    return QSize(500, 900);
  }

  void IDEvaluationBase::resetZoom()
  {
    spec_1d_->canvas()->resetZoom();
  }

  void IDEvaluationBase::setIntensityMode(int index)
  {
    intensity_button_group_->button(index)->setChecked(true);
    OpenMS::SpectrumCanvas::IntensityModes mode = (OpenMS::SpectrumCanvas::IntensityModes) index;
    switch (mode)
    {
    case OpenMS::SpectrumCanvas::IM_NONE:
    case OpenMS::SpectrumCanvas::IM_SNAP:
      spec_1d_->yAxis()->setLegend("# PSMs");
      break;

    case OpenMS::SpectrumCanvas::IM_LOG:
      spec_1d_->yAxis()->setLegend("# PSMs (log)");
      break;

    case OpenMS::SpectrumCanvas::IM_PERCENTAGE:
      spec_1d_->yAxis()->setLegend("PSMs [%]");
      break;
    }
    spec_1d_->setIntensityMode(mode);
  }

  bool IDEvaluationBase::getPoints(std::vector<PeptideIdentification>& peptides /* cannot be const, to avoid copy */,
                                   const std::vector<double>& q_value_thresholds, MSSpectrum& points)
  {
    points.clear(true);

    FalseDiscoveryRate fdr;
    fdr.setParameters(param_.copy("fdr:", true));
    try
    {
      fdr.apply(peptides); // computes a q-value (if its params are correct)
    }
    catch (Exception::MissingInformation)
    {
      LOG_FATAL_ERROR << "Tool failed due to missing information (see above)." << std::endl;
      return false;
    }

    // get list of q-values and sort them
    std::vector<double> q_values;
    q_values.reserve(peptides.size());
    for (vector<PeptideIdentification>::iterator it = peptides.begin(); it != peptides.end(); ++it)
    {
      it->assignRanks();
      if (it->getHits().size() > 0)
        q_values.push_back(it->getHits()[0].getScore());
    }
    std::sort(q_values.begin(), q_values.end());


    for (Size i = 0; i < q_value_thresholds.size(); ++i)
    {
      // get position in sorted q-values where cutoff is reached
      std::vector<double>::iterator pos = std::upper_bound(q_values.begin(), q_values.end(), q_value_thresholds[i]);
      Peak1D p;
      p.setMZ(q_value_thresholds[i] * 100);
      p.setIntensity(std::distance(q_values.begin(), pos));
      points.push_back(p);
    }

    return true;
  }

  void IDEvaluationBase::openFileDialog()
  {
    QString file_name = QFileDialog::getOpenFileName(this, tr("Open search result"), current_path_.toQString(), tr("search result (*.idXML)"));

    addSearchFile(file_name);
  }

  bool IDEvaluationBase::loadFiles(const StringList& list)
  {
    bool good = true;
    for (StringList::const_iterator it = list.begin(); it != list.end(); ++it)
    {
      if (!addSearchFile(*it)) good = false;
    }
    return good;
  }

  void IDEvaluationBase::setVisibleArea(double low, double high)
  {
    DRange<2> range(low * 100, -1, high * 100, -1);
    spec_1d_->canvas()->setVisibleArea(range);
  }


  bool IDEvaluationBase::loadCurve(const String& file_name, MSSpectrum& points)
  {
    if (FileHandler::getType(file_name) != FileTypes::IDXML)
    {
      LOG_ERROR << "The file '" << file_name << "' is not an .idXML file" << std::endl;
      return false;
    }

    std::vector<ProteinIdentification> prot_ids;
    std::vector<PeptideIdentification> pep_ids;
    IdXMLFile().load(file_name, prot_ids, pep_ids);
    String ln = pep_ids[0].getScoreType(); // grab name here, since FDR-calculation will overwrite it with "q-value"
    bool ret = getPoints(pep_ids, q_value_thresholds_, points); // FDR calculation failed?
    points.setMetaValue("search_engine", ln);

    return ret; 
  }

  bool IDEvaluationBase::addSearchFile(const String& file_name)
  {
    MSSpectrum points;
    if (!loadCurve(file_name, points)) return false;

    data_.addSpectrum(points);

    PeakMap* exp = new PeakMap();
    exp->addSpectrum(points);
    spec_1d_->canvas()->addLayer(SpectrumCanvas::ExperimentSharedPtrType(exp));
    spec_1d_->canvas()->setLayerName(spec_1d_->canvas()->getLayerCount() - 1, points.getMetaValue("search_engine"));
    // set intensity mode (after spectrum has been added!)
    setIntensityMode((int) SpectrumCanvas::IM_SNAP);

    return true;
  }
  
  const PeakMap& IDEvaluationBase::getPoints() const
  {
     return data_;
  }

  void IDEvaluationBase::saveImageAs()
  {
    QString cp = current_path_.toQString();
    QString file_name = QFileDialog::getSaveFileName(this, tr("Save image"), cp, tr("Images (*.svg *.png *.jpg)"));
    String error;
    if (!exportAsImage(file_name, error))
    {
      QMessageBox::warning(this, tr("Error"),
                           tr("Unable to save image to \n") +
                           file_name);
    }
  }

  StringList IDEvaluationBase::getSupportedImageFormats()
  {
    return ListUtils::create<String>("png,jpg,svg"); // make sure this is lower-case
  }

  bool IDEvaluationBase::exportAsImage(const QString& file_name, String& error_message, const QString& format)
  {
    if (file_name == "")
    {
      error_message = "Empty filename given!";
      return false;
    }

    QString suffix = format;
    if (suffix == "")
    {
      suffix = file_name.right(3);
    }

    if (!(suffix.compare("svg", Qt::CaseInsensitive) == 0) &&
        !(suffix.compare("png", Qt::CaseInsensitive) == 0) &&
        !(suffix.compare("jpg", Qt::CaseInsensitive) == 0))
    {
      error_message = "Unsupported format given('" + suffix + "')!";
      return false;
    }
    bool svg = (suffix.compare("svg", Qt::CaseInsensitive) == 0);

    // QSize items_bounding_rect = spec_1d_->size();
    // qreal wh_proportion = (qreal)(items_bounding_rect.width()) / (qreal)(items_bounding_rect.height());
    // bool w_larger_than_h = wh_proportion > 1;

    // qreal small_edge_length = svg ? 500 : 4000;

    // qreal x2, y2;
    // if (w_larger_than_h)
    // {
    //   x2 = wh_proportion * small_edge_length;
    //   y2 = small_edge_length;
    // }
    // else
    // {
    //   x2 = small_edge_length;
    //   y2 = (1.0 / wh_proportion) * small_edge_length;
    // }
    
    double h = param_.getValue("image:height");
    double w = param_.getValue("image:width");
    setGeometry(QRect(0, 0, w, h)); // does the layout -- otherwise we'd need show() to get it right

    if (svg)
    {
      QSvgGenerator svg_gen;
      svg_gen.setFileName(file_name);
      //svg_gen.setSize(items_bounding_rect.width(), items_bounding_rect.height());
      //svg_gen.setViewBox(QRect(x1, y1, x2, y2));
      svg_gen.setTitle(tr("Title (TBD)"));
      svg_gen.setDescription(tr("Description (TBD)"));
      QPainter painter(&svg_gen);
      spec_1d_->renderForImage(painter); //, QRectF(), items_bounding_rect);
    }
    else
    {
      spec_1d_->resize(w, h);

      QImage img(w, h, QImage::Format_ARGB32_Premultiplied);
      QPainter painter(&img);
      spec_1d_->renderForImage(painter); //, QRectF(), items_bounding_rect);
      painter.end();
      bool r = img.save(file_name, format.toStdString().c_str());
      if (!r)
      {
        error_message = "Could not save image to '" + file_name + "' with format '" + format + "'!";
        return false;
      }
    }
    return true;
  }

  void IDEvaluationBase::showURL()
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

  void IDEvaluationBase::showStatusMessage(string msg, OpenMS::UInt time)
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

  void IDEvaluationBase::showAboutDialog()
  {
    //dialog and grid layout
    QDialog* dlg = new QDialog(this);
    QGridLayout* grid = new QGridLayout(dlg);
    dlg->setWindowTitle("About IDEvaluation");

    QLabel* label = new QLabel(dlg);
    label->setPixmap(QPixmap(":/TOPP_about.png"));
    grid->addWidget(label, 0, 0);

    //text
    QString text = QString("<BR>"
                           "<FONT size=+3>IDEvaluation</font><BR>"
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
                             VersionInfo::getRevision() != "" ? QString(" (") + VersionInfo::getRevision().toQString() + ")" : "");    QLabel* text_label = new QLabel(text, dlg);
    grid->addWidget(text_label, 0, 1, Qt::AlignTop | Qt::AlignLeft);

    //execute
    dlg->exec();
  }

  void IDEvaluationBase::showLogMessage_(IDEvaluationBase::LogState state, const String& heading, const String& body)
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

  void IDEvaluationBase::keyPressEvent(QKeyEvent* e)
  {
    if (e->key() == Qt::Key_F5)
    {
      /*TOPPASWidget* tw = activeWindow_();
      if (!tw)
      {
          e->ignore();
          return;
      }
      TOPPASScene* ts = tw->getScene();
      ts->runPipeline();
      e->accept();
*/
    }
    e->ignore();
  }

  void IDEvaluationBase::closeEvent(QCloseEvent* /* event */)
  {
    //ws_->closeAllWindows();
    //event->accept();
  }

} //namespace OpenMS

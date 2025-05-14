// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/Plot1DWidget.h>
#include <OpenMS/VISUAL/Plot2DWidget.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data
#include <QApplication>

#include <boost/make_shared.hpp>

using namespace OpenMS;
using namespace std;

/**
 * @brief Loads a single mass spectrum from a file and displays it in a 1D plot using OpenMS and Qt.
 *
 * Initializes a Qt application, loads a spectrum from a DTA file, adds it to an annotated experiment,
 * and visualizes the spectrum in a 1D plot widget. The application window is displayed and the Qt event loop is started.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return The exit code from the Qt application event loop.
 */
Int main(int argc, const char** argv)
{
  String tutorial_data_path(OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_Spectrum1D.dta"));

  QApplication app(argc, const_cast<char**>(argv));

  AnnotatedMSRun exp;
  auto exp_sptr = boost::make_shared<AnnotatedMSRun>();
  MSSpectrum spec;
  // demonstrating how to load a single spectrum from file formats which only contain a single spec
  // alternatively: use FileHandler().loadExperiment() if you need an experiment anyway
  FileHandler().loadSpectrum(tutorial_data_path, spec, {FileTypes::DTA});
  exp_sptr->getMSExperiment().addSpectrum(spec);
  LayerDataBase::ODExperimentSharedPtrType on_disc_exp_sptr(new OnDiscMSExperiment());
  Plot1DWidget widget(Param(), DIM::Y, nullptr);
  widget.canvas()->addPeakLayer(exp_sptr, on_disc_exp_sptr);
  widget.show();

  return app.exec();
} // end of main

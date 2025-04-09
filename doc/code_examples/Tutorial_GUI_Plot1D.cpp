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

using namespace OpenMS;
using namespace std;

Int main(int argc, const char** argv)
{
  String tutorial_data_path(OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_Spectrum1D.dta"));

  QApplication app(argc, const_cast<char**>(argv));

  PeakMap exp;
  MSSpectrum spec;
  // demonstrating how to load a single spectrum from file formats which only contain a single spec
  // alternatively: use FileHandler().loadExperiment() if you need an experiment anyway
  FileHandler().loadSpectrum(tutorial_data_path, spec, {FileTypes::DTA});
  exp.addSpectrum(spec);
  LayerDataBase::ExperimentSharedPtrType exp_sptr(new PeakMap(exp));
  LayerDataBase::ODExperimentSharedPtrType on_disc_exp_sptr(new OnDiscMSExperiment());
  Plot1DWidget widget(Param(), DIM::Y, nullptr);
  widget.canvas()->addPeakLayer(exp_sptr, on_disc_exp_sptr);
  widget.show();

  return app.exec();
} // end of main

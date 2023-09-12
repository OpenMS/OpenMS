// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/Plot1DWidget.h>
#include <OpenMS/VISUAL/Plot2DWidget.h>
#include <OpenMS/openms_data_path.h>
#include <QApplication>

using namespace OpenMS;
using namespace std;

Int main(int argc, const char** argv)
{
  String tutorial_data_path(OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_Spectrum1D.dta"));
  if (argc >= 2)
  { // the path to the data can be given on the command line
    tutorial_data_path = argv[1];
  }

  QApplication app(argc, const_cast<char**>(argv));

  PeakMap exp;
  exp.resize(1);
  DTAFile().load(tutorial_data_path, exp[0]);
  LayerDataBase::ExperimentSharedPtrType exp_sptr(new PeakMap(exp));
  LayerDataBase::ODExperimentSharedPtrType on_disc_exp_sptr(new OnDiscMSExperiment());
  auto* widget = new Plot1DWidget(Param(), DIM::Y, nullptr);
  widget->canvas()->addPeakLayer(exp_sptr, on_disc_exp_sptr);
  widget->show();

  return app.exec();
} // end of main

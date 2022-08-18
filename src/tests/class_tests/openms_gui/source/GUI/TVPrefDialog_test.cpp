// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <TVPrefDialog_test.h>

#include <qcheckbox.h>

using namespace OpenMS;

// delay in ms
// higher values together with 'dialog_.show();' can be useful for debugging this test
constexpr int DELAY {15};

void testCheckBox(QCheckBox* cb)
{
  qInfo() << qPrintable("Testing Checkbox: " + cb->objectName());
  Qt::CheckState prev_state = cb->checkState();
  QTest::mouseClick(cb, Qt::LeftButton, 0, cb->rect().bottomLeft());
  QTest::qWait(DELAY);
  QVERIFY2(prev_state != cb->checkState(), qPrintable("'" + cb->objectName() + "' check box didn't change its state by clicking on it."));
}

void TestTVPrefDialog::testConstruction()
{
  QVERIFY2(UI->tabWidget, "Tab widget not created.");
  QVERIFY2(UI->buttonBox, "Dialog button box not created.");

  // General tab
  QVERIFY2(UI->browse_default, "'Browse' path button not created.");
  QVERIFY2(UI->default_path, "Default path line edit not created.");
  QVERIFY2(UI->browse_plugins, "'Browse' plugins button not created.");
  QVERIFY2(UI->plugins_path, "Plugins path line edit not created.");
  QVERIFY2(UI->default_path_current, "'Use current path' check box not created.");
  QVERIFY2(UI->map_default, "'Default map view' combo box not created.");
  QVERIFY2(UI->map_cutoff, "'Low intensity cutoff' combo box not created.");
  QVERIFY2(UI->on_file_change, "'Action when file changes' combo box not created.");
  QVERIFY2(UI->use_cached_ms1, "'Cache MS1 spectra to disk' check box not created.");
  QVERIFY2(UI->use_cached_ms2, "'Cache MS2 spectra to disk' check box not created.");
  QVERIFY2(UI->default_path_label, "'Default path:' label not created.");
  QVERIFY2(UI->plugins_path_label, "'Plugins path:' label not created.");
  QVERIFY2(UI->default_map_label, "'Default map view:' label not created.");
  QVERIFY2(UI->intensity_cutoff_label, "'Low intensity cutoff:' label not created.");
  QVERIFY2(UI->file_change_label, "'Action when file changes:' label not created.");
  QVERIFY2(UI->caching_label, "'Caching Strategy:' label not created.");

  // 1D tab
  QVERIFY2(UI->color_1D, "Peak color color selector not created.");
  QVERIFY2(UI->selected_1D, "Selected peak color color selector not created.");
  QVERIFY2(UI->icon_1D, "Icon color color selector not created.");
  QVERIFY2(UI->peak_1D_label, "'Peak color:' label not created.");
  QVERIFY2(UI->selected_1D_label, "'Selected peak color:' label not created.");
  QVERIFY2(UI->icon_1D_label, "'Icon color:' label not created.");

  // 2D tab
  QVERIFY2(UI->peak_2D, "Peak 2D multi gradient selector not created.");
  QVERIFY2(UI->mapping_2D, "'m/z axis' combo box not created.");
  QVERIFY2(UI->feature_icon_2D, "'feature icon' combo box not created.");
  QVERIFY2(UI->feature_icon_size_2D, "'feature icon size' spin box not created.");
  QVERIFY2(UI->peak_2D_label, "'Peak gradient:' label not created.");
  QVERIFY2(UI->mz_axis_label, "'m/z axis:' label not created.");
  QVERIFY2(UI->feature_icon_2D_label, "'feature icon:' label not created.");
  QVERIFY2(UI->feature_icon_size_2D_label, "'feature icon size:' label not created.");

  // 3D tab
  QVERIFY2(UI->peak_3D, "Peak 3D multi gradient selector not created.");
  QVERIFY2(UI->shade_3D, "Shade mode combo box not created.");
  QVERIFY2(UI->line_width_3D, "Line width spin box not created.");
  QVERIFY2(UI->peak_3D_label, "'Peak colors:' label not created.");
  QVERIFY2(UI->shade_3D_label, "'Shade mode:' label not created.");
  QVERIFY2(UI->line_width_3D_label, "'Line width:' label not created.");

  // TSG tab
  QVERIFY2(UI->param_editor_spec_gen_, "Parameter editor not created.");
  QVERIFY2(UI->tolerance, "Tolerance double spin box not created.");
  QVERIFY2(UI->unit, "Unit combo box not created.");
  QVERIFY2(UI->tolerance_label, "'Alignment Tolerance:' label not created.");
  QVERIFY2(UI->param_usage_label, "'These settings aare used to ...' label not created.");
}

void TestTVPrefDialog::testGui()
{
  // TODO: currently Qt::ComboBox, OpenMS::MultiGradientSelector and OpenMS::ParamEditor are not tested!
  // I was unable to get the positions of the interactable parts of these objects for a 'QTest::mouseClick'

  enum
  {
    GENERAL,
    ONE_D,
    TWO_D,
    THREE_D,
    TSG
  };

  QTabBar* tab_bar = UI->tabWidget->tabBar();

  /////////////////////////////////////////////
  //             'General' tab               //
  /////////////////////////////////////////////

  // click on tab
  QTest::mouseClick(tab_bar, Qt::LeftButton, 0, tab_bar->tabRect(GENERAL).center());
  QTest::qWait(DELAY);

  // test path input methods
  
  // file dialog (default path)
  QTimer::singleShot(DELAY, this, &TestTVPrefDialog::checkFileDialog_);
  QTest::mouseClick(UI->browse_default, Qt::LeftButton, 0, QPoint());
  QTest::qWait(DELAY);

  // line edit (default path)
  UI->default_path->clear();
  QTest::keyClicks(UI->default_path, "C:\\dev");
  QVERIFY2(UI->default_path->text() == "C:\\dev", "Line edit for default path broken.");
  QTest::qWait(DELAY);
  UI->default_path->clear();

  // file dialog (plugins path)
  QTimer::singleShot(DELAY, this, &TestTVPrefDialog::checkFileDialog_);
  QTest::mouseClick(UI->browse_plugins, Qt::LeftButton, 0, QPoint());
  QTest::qWait(DELAY);

  // line edit (plugins path)
  UI->plugins_path->clear();
  QTest::keyClicks(UI->plugins_path, "C:\\dev");
  QVERIFY2(UI->plugins_path->text() == "C:\\dev", "Line edit for plugins path broken.");
  QTest::qWait(DELAY);
  UI->plugins_path->clear();


  // test check boxes
  testCheckBox(UI->default_path_current);
  testCheckBox(UI->use_cached_ms1);
  testCheckBox(UI->use_cached_ms2);

  // TODO: 'map_default' (Qt::ComboBox), 'map_cutoff' (Qt::ComboBox), 'on_file_change' (Qt::ComboBox)

  /////////////////////////////////////////////
  //             '1D view' tab               //
  /////////////////////////////////////////////

  // click on tab
  QTest::mouseClick(tab_bar, Qt::LeftButton, 0, tab_bar->tabRect(ONE_D).center());
  QTest::qWait(DELAY);

  // test color selectors
  // Opening a dialog stops the test until it's closed. Therefore use 'singleShot' to call a function with a delay.
  QTimer::singleShot(DELAY, this, &TestTVPrefDialog::checkColorDialog_);
  QTest::mouseClick(UI->color_1D, Qt::LeftButton, 0, QPoint());
  QTest::qWait(DELAY);

  QTimer::singleShot(DELAY, this, &TestTVPrefDialog::checkColorDialog_);
  QTest::mouseClick(UI->selected_1D, Qt::LeftButton, 0, QPoint());
  QTest::qWait(DELAY);

  QTimer::singleShot(DELAY, this, &TestTVPrefDialog::checkColorDialog_);
  QTest::mouseClick(UI->icon_1D, Qt::LeftButton, 0, QPoint());
  QTest::qWait(DELAY);

  /////////////////////////////////////////////
  //             '2D view' tab               //
  /////////////////////////////////////////////

  // click on tab
  QTest::mouseClick(tab_bar, Qt::LeftButton, 0, tab_bar->tabRect(TWO_D).center());
  QTest::qWait(DELAY);

  // test spin box
  UI->feature_icon_size_2D->clear();
  QTest::keyClicks(UI->feature_icon_size_2D, "5");
  QTest::qWait(DELAY);
  QVERIFY(5 == UI->feature_icon_size_2D->value());

  // TODO:: 'peak_2D' (OpenMS::MultiGradientSelector), 'mapping_2D' (Qt::ComboBox), 'feature_icon_2D' (Qt::ComboBox)
  
  /////////////////////////////////////////////
  //             '3D view' tab               //
  /////////////////////////////////////////////

  // click on tab
  QTest::mouseClick(tab_bar, Qt::LeftButton, 0, tab_bar->tabRect(THREE_D).center());
  QTest::qWait(DELAY);

  // test spin box
  UI->line_width_3D->clear();
  QTest::keyClicks(UI->line_width_3D, "2");
  QTest::qWait(DELAY);
  QVERIFY(2 == UI->line_width_3D->value());

  // TODO: 'shade_3D' (Qt::ComboBox), 'peak_3D' (OpenMS::MultiGradientSelector)

  /////////////////////////////////////////////
  //               'TSG' tab                 //
  /////////////////////////////////////////////

  // click on tab
  QTest::mouseClick(tab_bar, Qt::LeftButton, 0, tab_bar->tabRect(TSG).center());
  QTest::qWait(DELAY);

  // test spin box
  UI->tolerance->clear();
  QTest::keyClicks(UI->tolerance, QLocale::system().toString(0.5));
  QTest::qWait(DELAY);
  QVERIFY(0.5 == UI->tolerance->value());

  // TODO: 'param_editor_spec_gen_' (OpenMS::ParamEditor), 'unit' (Qt::ComboBox)
}

void OpenMS::TestTVPrefDialog::testParamExport()
{
  // check default parameters
  Param dialog_param = dialog_.getParam(); 
  
#define PARAM(a) dialog_param.getValue(a)

  QVERIFY2(PARAM("default_path") == "", "'Default path' param value not exported correctly.");
  QVERIFY2(PARAM("default_path_current") == "false", "'Use current path' param value not exported correctly.");
  QVERIFY2(PARAM("use_cached_ms1") == "false", "'Cache ms1 spectra' param value not exported correctly.");
  QVERIFY2(PARAM("use_cached_ms2") == "false", "'Cache ms2 spectra' param value not exported correctly.");
  QVERIFY2(PARAM("default_map_view") == "2d", "'Default map view' param value not exported correctly.");
  QVERIFY2(PARAM("intensity_cutoff") == "on", "'Itensity cutoff' param value not exported correctly.");
  QVERIFY2(PARAM("on_file_change") == "none", "'Action when file changes' param value not exported correctly.");

  QVERIFY2(PARAM("1d:peak_color") == "#ffffff", "'1D peak color' param value not exported correctly.");
  QVERIFY2(PARAM("1d:highlighted_peak_color") == "#ffffff", "'1D highlighted peak color' param value not exported correctly.");
  QVERIFY2(PARAM("1d:icon_color") == "#ffffff", "'1D icon color' param value not exported correctly.");

  QVERIFY2(PARAM("2d:dot:gradient") == "Linear|0,#ffffff;100,#000000", "'2D peak gradient' param value not exported correctly.");
  QVERIFY2(PARAM("2d:mapping_of_mz_to") == "x_axis", "'2D mapping of mz to' param value not exported correctly.");
  QVERIFY2(PARAM("2d:dot:feature_icon") == "diamond", "'2D feature icon' param value not exported correctly.");
  QVERIFY2(int(PARAM("2d:dot:feature_icon_size")) == 3, "'2D feature icon size' param value not exported correctly.");

  QVERIFY2(PARAM("3d:dot:gradient") == "Linear|0,#ffffff;100,#000000", "'3D peak gradient' param value not exported correctly.");
  QVERIFY2(int(PARAM("3d:dot:shade_mode")) == 0, "'3D shade mode' param value not exported correctly.");
  QVERIFY2(int(PARAM("3d:dot:line_width")) == 1, "'3D line width' param value not exported correctly.");

  QVERIFY2(PARAM("idview:tsg:isotope_model") == "none", "TSG: 'isotope model' param value not exported correctly.");
  QVERIFY2(int(PARAM("idview:tsg:max_isotope")) == 2, "TSG: 'max isotope' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:max_isotope_probability")) == 0.05, "TSG: 'max isotope probability' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_metainfo") == "false", "TSG: 'add metainfo' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_losses") == "false", "TSG: 'add losses' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:sort_by_position") == "true", "TSG: 'sort by position' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_precursor_peaks") == "false", "TSG: 'add precursor peaks' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_all_precursor_charges") == "false", "TSG: 'add all precursor charges' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_abundant_immonium_ions") == "false", "TSG: 'add abundant immonium ions' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_first_prefix_ion") == "false", "TSG: 'add first prefix ion' param value not exported correctly.");
  for (String s : {'a', 'c', 'x', 'x', 'z'})
  {
    QVERIFY2(PARAM("idview:tsg:add_" + s + "_ions") == "false", qPrintable("TSG: 'add " + s.toQString() + " ions' param value not exported correctly."));
    QVERIFY2(double(PARAM("idview:tsg:" + s + "_intensity")) == 1.0, qPrintable("TSG: '" + s.toQString() + " intensity' param value not exported correctly."));
  }
  QVERIFY2(PARAM("idview:tsg:add_b_ions") == "true", "TSG: 'add b ions' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:b_intensity")) == 1.0, "TSG: 'b intensity' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_y_ions") == "true", "TSG: 'add y ions' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:y_intensity")) == 1.0, "TSG: 'y intensity' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:relative_loss_intensity")) == 0.1, "TSG: 'relative loss intensity' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:precursor_intensity")) == 1.0, "TSG: 'precursor intensity' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:precursor_H2O_intensity")) == 1.0, "TSG: 'precursor H2O intensity' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:precursor_NH3_intensity")) == 1.0, "TSG: 'precursor NH3 intensity' param value not exported correctly.");
  
  // set some custom parameters (different from the default)
  
  // General
  UI->default_path->setText("C:\\dev");                          // 'C:\dev' default path
  UI->default_path_current->setCheckState(Qt::CheckState::Checked); // use current path
  UI->map_default->setCurrentIndex(1);                          // 3D
  UI->map_cutoff->setCurrentIndex(1);                           // cut off 'off'
  UI->on_file_change->setCurrentIndex(2);                       // 'update automtically'
  UI->use_cached_ms1->setCheckState(Qt::CheckState::Checked);   // cache ms1
  UI->use_cached_ms2->setCheckState(Qt::CheckState::Checked);   // cache ms2

  // 1D
  UI->color_1D->setColor(QColor("#ff0000"));        // peak 1D 'red' ("#ff0000")
  UI->selected_1D->setColor(QColor("#8b0000"));     // selected peak 1D 'darkred' ("#8b0000")
  UI->icon_1D->setColor(QColor("#660000"));         // icon 1D 'crimson' ("#660000")

  // 2D
  UI->peak_2D->gradient().fromString("Linear|0,#ffaa00;6,#ff0000;14,#aa00ff;23,#5500ff;100,#000000"); // peak 2D: orange - red - purple - blue - black
  UI->mapping_2D->setCurrentIndex(1);       // y-axis
  UI->feature_icon_2D->setCurrentIndex(2);  // circle
  UI->feature_icon_size_2D->setValue(5);    // size: 5

  // 3D
  UI->peak_3D->gradient().fromString("Linear|0,#ffea00;6,#ffaa00;100,#ff0000"); // peak 3D: yellow - orange - red
  UI->shade_3D->setCurrentIndex(1); // smooth (index 1)
  UI->line_width_3D->setValue(2);   // line width: 2

  // TSG
  Param p;
  p.setValue("isotope_model", "coarse");
  p.setValue("max_isotope", 1);
  p.setValue("max_isotope_probability", 0.01);
  p.setValue("add_metainfo", "true");
  p.setValue("add_losses", "true");
  p.setValue("sort_by_position", "false");
  p.setValue("add_precursor_peaks", "true");
  p.setValue("add_all_precursor_charges", "true");
  p.setValue("add_abundant_immonium_ions", "true");
  p.setValue("add_first_prefix_ion", "true");
  for (String s : {'a', 'c', 'x', 'x', 'z'})
  {
    p.setValue("add_" + s + "_ions", "true");
    p.setValue(s + "_intensity", 0.8);
  }
  p.setValue("add_b_ions", "false");
  p.setValue("b_intensity", 0.8);
  p.setValue("add_y_ions", "false");
  p.setValue("y_intensity", 0.8);
  p.setValue("relative_loss_intensity", 0.2);
  p.setValue("precursor_intensity", 0.99);
  p.setValue("precursor_H2O_intensity", 0.95);
  p.setValue("precursor_NH3_intensity", 0.9);

  dialog_.setParam(p);
  UI->tolerance->setValue(0.2); // tolerance: 0.2
  UI->unit->setCurrentIndex(0); // Dalton

  // accept the dialog to save the parameters
  QTest::mouseClick(UI->buttonBox->button(QDialogButtonBox::Ok), Qt::LeftButton);

  // get parameters
  dialog_param = dialog_.getParam();

  // check validity
  QVERIFY2(PARAM("default_path") == "C:\\dev", "'Default path' param value not exported correctly.");
  QVERIFY2(PARAM("default_path_current") == "true", "'Use current path' param value not exported correctly.");
  QVERIFY2(PARAM("use_cached_ms1") == "true", "'Cache ms1 spectra' param value not exported correctly.");
  QVERIFY2(PARAM("use_cached_ms2") == "true", "'Cache ms2 spectra' param value not exported correctly.");
  QVERIFY2(PARAM("default_map_view") == "3d", "'Default map view' param value not exported correctly.");
  QVERIFY2(PARAM("intensity_cutoff") == "off", "'Itensity cutoff' param value not exported correctly.");
  QVERIFY2(PARAM("on_file_change") == "update automatically", "'Action when file changes' param value not exported correctly.");

  QVERIFY2(PARAM("1d:peak_color") == "#ff0000", "'1D peak color' param value not exported correctly.");
  QVERIFY2(PARAM("1d:highlighted_peak_color") == "#8b0000", "'1D highlighted peak color' param value not exported correctly.");
  QVERIFY2(PARAM("1d:icon_color") == "#660000", "'1D icon color' param value not exported correctly.");

  QVERIFY2(PARAM("2d:dot:gradient") == "Linear|0,#ffaa00;6,#ff0000;14,#aa00ff;23,#5500ff;100,#000000", "'2D peak gradient' param value not exported correctly.");
  QVERIFY2(PARAM("2d:mapping_of_mz_to") == "y_axis", "'2D mapping of mz to' param value not exported correctly.");
  QVERIFY2(PARAM("2d:dot:feature_icon") == "circle", "'2D feature icon' param value not exported correctly.");
  QVERIFY2(int(PARAM("2d:dot:feature_icon_size")) == 5, "'2D feature icon size' param value not exported correctly.");

  QVERIFY2(PARAM("3d:dot:gradient") == "Linear|0,#ffea00;6,#ffaa00;100,#ff0000", "'3D peak gradient' param value not exported correctly.");
  QVERIFY2(int(PARAM("3d:dot:shade_mode")) == 1, "'3D shade mode' param value not exported correctly.");
  QVERIFY2(int(PARAM("3d:dot:line_width")) == 2, "'3D line width' param value not exported correctly.");

  QVERIFY2(PARAM("idview:tsg:isotope_model") == "coarse", "TSG: 'isotope model' param value not exported correctly.");
  QVERIFY2(int(PARAM("idview:tsg:max_isotope")) == 1, "TSG: 'max isotope' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:max_isotope_probability")) == 0.01, "TSG: 'max isotope probability' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_metainfo") == "true", "TSG: 'add metainfo' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_losses") == "true", "TSG: 'add losses' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:sort_by_position") == "false", "TSG: 'sort by position' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_precursor_peaks") == "true", "TSG: 'add precursor peaks' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_all_precursor_charges") == "true", "TSG: 'add all precursor charges' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_abundant_immonium_ions") == "true", "TSG: 'add abundant immonium ions' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_first_prefix_ion") == "true", "TSG: 'add first prefix ion' param value not exported correctly.");
  for (String s : {'a', 'c', 'x', 'x', 'z'})
  {
    QVERIFY2(PARAM("idview:tsg:add_" + s + "_ions") == "true", qPrintable("TSG: 'add " + s.toQString() + " ions' param value not exported correctly."));
    QVERIFY2(double(PARAM("idview:tsg:" + s + "_intensity")) == 0.8, qPrintable("TSG: '" + s.toQString() + " intensity' param value not exported correctly."));
  }
  QVERIFY2(PARAM("idview:tsg:add_b_ions") == "false", "TSG: 'add b ions' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:b_intensity")) == 0.8, "TSG: 'b intensity' param value not exported correctly.");
  QVERIFY2(PARAM("idview:tsg:add_y_ions") == "false", "TSG: 'add y ions' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:y_intensity")) == 0.8, "TSG: 'y intensity' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:relative_loss_intensity")) == 0.2, "TSG: 'relative loss intensity' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:precursor_intensity")) == 0.99, "TSG: 'precursor intensity' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:precursor_H2O_intensity")) == 0.95, "TSG: 'precursor H2O intensity' param value not exported correctly.");
  QVERIFY2(double(PARAM("idview:tsg:precursor_NH3_intensity")) == 0.9, "TSG: 'precursor NH3 intensity' param value not exported correctly.");
}

void OpenMS::TestTVPrefDialog::checkFileDialog_()
{
  // get the active window
  QWidget* active_widget = QApplication::activeModalWidget();
  if (active_widget->inherits("QFileDialog")) // if it's a file dialog, close it
  {
    QFileDialog* fd = qobject_cast<QFileDialog*>(active_widget);
    fd->close(); // for some reason closing it with 'Qt::Key_Enter' doesn't work
    qInfo() << "Closing File Dialog.";
    QVERIFY(true);
    return;
  }
  QVERIFY(false);
}
void OpenMS::TestTVPrefDialog::checkColorDialog_()
{
  // get the active window
  QWidget* active_widget = QApplication::activeModalWidget();
  if (active_widget->inherits("QColorDialog")) // if it's a color dialog, close it
  {
    QColorDialog* cd = qobject_cast<QColorDialog*>(active_widget);
    QTest::keyClick(cd, Qt::Key_Enter);
    qInfo() << "Closing Color Dialog.";
    QVERIFY(true);
    return;
  }
  QVERIFY(false);
}

// expands to a simple main() method that runs all the private slots
QTEST_MAIN(TestTVPrefDialog)
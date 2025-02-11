// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPViewPrefDialog.h>
#include <ui_TOPPViewPrefDialog.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <QtWidgets/QFileDialog>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    TOPPViewPrefDialog::TOPPViewPrefDialog(QWidget* parent) :
      QDialog(parent),
      ui_(new Ui::TOPPViewPrefDialogTemplate),
      tsg_param_(TheoreticalSpectrumGenerator().getParameters())
    {
      ui_->setupUi(this);
      ui_->param_editor_spec_gen_->load(tsg_param_);
      connect(ui_->browse_default, &QPushButton::clicked, this, &TOPPViewPrefDialog::browseDefaultPath_);
      connect(ui_->browse_plugins, &QPushButton::clicked, this, &TOPPViewPrefDialog::browsePluginsPath_);
    }

    TOPPViewPrefDialog::~TOPPViewPrefDialog()
    {
      delete ui_;
    }

    const char* tsg_prefix = "idview:tsg:";

    void TOPPViewPrefDialog::setParam(const Param& param)
    {
      param_ = getParam(); // get our own defaults

      // make sure the params we write (using getParam) are the same as the ones we can read
      param_.update(param, true, true, true, true, OpenMS_Log_info);

      // general tab
      ui_->default_path->setText(String(param_.getValue("default_path").toString()).toQString());
      ui_->default_path_current->setChecked(param_.getValue("default_path_current").toBool());
      ui_->plugins_path->setText(String(param_.getValue("plugins_path").toString()).toQString());
      ui_->use_cached_ms1->setChecked(param_.getValue("use_cached_ms1").toBool());
      ui_->use_cached_ms2->setChecked(param_.getValue("use_cached_ms2").toBool());

      ui_->map_default->setCurrentIndex(ui_->map_default->findText(String(param_.getValue("default_map_view").toString()).toQString()));
      ui_->map_cutoff->setCurrentIndex(ui_->map_cutoff->findText(String(param_.getValue("intensity_cutoff").toString()).toQString()));
      ui_->on_file_change->setCurrentIndex(ui_->on_file_change->findText(String(param_.getValue("on_file_change").toString()).toQString()));

      // 1D view
      ui_->color_1D->setColor(QColor(String(param_.getValue("1d:peak_color").toString()).toQString()));
      ui_->selected_1D->setColor(QColor(String(param_.getValue("1d:highlighted_peak_color").toString()).toQString()));
      ui_->icon_1D->setColor(QColor(String(param_.getValue("1d:icon_color").toString()).toQString()));

      // 2D view
      ui_->peak_2D->gradient().fromString(param_.getValue("2d:dot:gradient"));
      ui_->feature_icon_2D->setCurrentIndex(ui_->feature_icon_2D->findText(String(param_.getValue("2d:dot:feature_icon").toString()).toQString()));
      ui_->feature_icon_size_2D->setValue((Int)param_.getValue("2d:dot:feature_icon_size"));

      // 3D view
      ui_->peak_3D->gradient().fromString(param_.getValue("3d:dot:gradient"));
      ui_->shade_3D->setCurrentIndex((Int)param_.getValue("3d:dot:shade_mode"));
      ui_->line_width_3D->setValue((Int)param_.getValue("3d:dot:line_width"));

      // TSG view
      tsg_param_ = param_.copy(tsg_prefix, true);
      ui_->param_editor_spec_gen_->load(tsg_param_);
      ui_->tolerance->setValue((double)param_.getValue("idview:align:tolerance"));
      ui_->unit->setCurrentIndex(ui_->unit->findText(String(param_.getValue("idview:align:is_relative_tolerance") == "true" ? "ppm" : "Da").toQString()));
    }

    String fromCheckState(const Qt::CheckState cs)
    {
      switch (cs)
      {
        case Qt::CheckState::Checked: return "true";
        case Qt::CheckState::Unchecked: return "false";
        default: throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Checkbox had unexpected state", String(cs));
      }
    }

    Param TOPPViewPrefDialog::getParam() const
    {
      Param p; 
      p.setValue("default_path", ui_->default_path->text().toStdString());
      p.setValue("default_path_current", fromCheckState(ui_->default_path_current->checkState()));

      p.setValue("plugins_path", ui_->plugins_path->text().toStdString());

      p.setValue("use_cached_ms1", fromCheckState(ui_->use_cached_ms1->checkState()));
      p.setValue("use_cached_ms2", fromCheckState(ui_->use_cached_ms2->checkState()));

      p.setValue("default_map_view", ui_->map_default->currentText().toStdString());
      p.setValue("intensity_cutoff", ui_->map_cutoff->currentText().toStdString());
      p.setValue("on_file_change", ui_->on_file_change->currentText().toStdString());

      p.setValue("1d:peak_color", ui_->color_1D->getColor().name().toStdString());
      p.setValue("1d:highlighted_peak_color", ui_->selected_1D->getColor().name().toStdString());
      p.setValue("1d:icon_color", ui_->icon_1D->getColor().name().toStdString());

      p.setValue("2d:dot:gradient", ui_->peak_2D->gradient().toString());
      p.setValue("2d:dot:feature_icon", ui_->feature_icon_2D->currentText().toStdString());
      p.setValue("2d:dot:feature_icon_size", ui_->feature_icon_size_2D->value());

      p.setValue("3d:dot:gradient", ui_->peak_3D->gradient().toString());
      p.setValue("3d:dot:shade_mode", ui_->shade_3D->currentIndex());
      p.setValue("3d:dot:line_width", ui_->line_width_3D->value());

      // TSG view
      ui_->param_editor_spec_gen_->store(); // to tsg_param_
      p.insert(tsg_prefix, tsg_param_);
      // from SpectrumAlignment
      p.setValue("idview:align:tolerance", ui_->tolerance->value(), "Alignment tolerance value");
      p.setValue("idview:align:is_relative_tolerance", ui_->unit->currentText().toStdString() == "ppm" ? "true" : "false", "Alignment tolerance unit (Da, ppm)");

      param_ = p;
        
      return param_;
    }

    void TOPPViewPrefDialog::browseDefaultPath_()
    {
      QString path = QFileDialog::getExistingDirectory(this, "Choose a directory", ui_->default_path->text());
      if (!path.isEmpty())
      {
        ui_->default_path->setText(path);
      }
    }

    void TOPPViewPrefDialog::browsePluginsPath_()
    {
      QString path = QFileDialog::getExistingDirectory(this, "Choose a directory", ui_->plugins_path->text());
      if (!path.isEmpty())
      {
        ui_->plugins_path->setText(path);
      }
    }

  }   //namespace Internal
} //namespace OpenMS

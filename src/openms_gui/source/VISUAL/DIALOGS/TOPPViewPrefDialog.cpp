// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPViewPrefDialog.h>
#include <ui_TOPPViewPrefDialog.h>


#include <QtWidgets/QFileDialog>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    TOPPViewPrefDialog::TOPPViewPrefDialog(QWidget* parent) :
      QDialog(parent),
      ui_(new Ui::TOPPViewPrefDialogTemplate)
    {
      ui_->setupUi(this);
      connect(ui_->browse_default, &QPushButton::clicked, this, &TOPPViewPrefDialog::browseDefaultPath_);
      connect(ui_->browse_temp, &QPushButton::clicked, this, &TOPPViewPrefDialog::browseTempPath_);
    }

    TOPPViewPrefDialog::~TOPPViewPrefDialog()
    {
      delete ui_;
    }

    void TOPPViewPrefDialog::setParam(const Param& param)
    {
      param_ = param;

      // --------------------------------------------------------------------
      // Set dialog entries from current parameter object (default values)

      // default
      ui_->default_path->setText(param_.getValue("preferences:default_path").toQString());
      ui_->default_path_current->setChecked(param_.getValue("preferences:default_path_current").toBool());
      ui_->use_cached_ms1->setChecked(param_.getValue("preferences:use_cached_ms1").toBool());
      ui_->use_cached_ms2->setChecked(param_.getValue("preferences:use_cached_ms2").toBool());
   
      ui_->temp_path->setText(param_.getValue("preferences:tmp_file_path").toQString());
      ui_->recent_files->setValue((Int)param_.getValue("preferences:number_of_recent_files"));
      ui_->map_default->setCurrentIndex(ui_->map_default->findText(param_.getValue("preferences:default_map_view").toQString()));
      ui_->map_cutoff->setCurrentIndex(ui_->map_cutoff->findText(param_.getValue("preferences:intensity_cutoff").toQString()));
      ui_->on_file_change->setCurrentIndex(ui_->on_file_change->findText(param_.getValue("preferences:on_file_change").toQString()));

      // 1D view
      ui_->color_1D->setColor(QColor(param_.getValue("preferences:1d:peak_color").toQString()));
      ui_->selected_1D->setColor(QColor(param_.getValue("preferences:1d:highlighted_peak_color").toQString()));
      ui_->icon_1D->setColor(QColor(param_.getValue("preferences:1d:icon_color").toQString()));

      // 2D view
      ui_->peak_2D->gradient().fromString(param_.getValue("preferences:2d:dot:gradient"));
      ui_->mapping_2D->setCurrentIndex(ui_->mapping_2D->findText(param_.getValue("preferences:2d:mapping_of_mz_to").toQString()));
      ui_->feature_icon_2D->setCurrentIndex(ui_->feature_icon_2D->findText(param_.getValue("preferences:2d:dot:feature_icon").toQString()));
      ui_->feature_icon_size_2D->setValue((Int)param_.getValue("preferences:2d:dot:feature_icon_size"));

      // 3D view
      ui_->peak_3D->gradient().fromString(param_.getValue("preferences:3d:dot:gradient"));
      ui_->shade_3D->setCurrentIndex((Int)param_.getValue("preferences:3d:dot:shade_mode"));
      ui_->line_width_3D->setValue((Int)param_.getValue("preferences:3d:dot:line_width"));

      // id view
      ui_->a_intensity->setValue((double)param_.getValue("preferences:idview:a_intensity"));
      ui_->b_intensity->setValue((double)param_.getValue("preferences:idview:b_intensity"));
      ui_->c_intensity->setValue((double)param_.getValue("preferences:idview:c_intensity"));
      ui_->x_intensity->setValue((double)param_.getValue("preferences:idview:x_intensity"));
      ui_->y_intensity->setValue((double)param_.getValue("preferences:idview:y_intensity"));
      ui_->z_intensity->setValue((double)param_.getValue("preferences:idview:z_intensity"));
      ui_->tolerance->setValue((double)param_.getValue("preferences:idview:tolerance"));

      ui_->relative_loss_intensity->setValue((double)param_.getValue("preferences:idview:relative_loss_intensity"));

      QList<QListWidgetItem*> a_ions = ui_->ions_list_widget->findItems("A-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> b_ions = ui_->ions_list_widget->findItems("B-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> c_ions = ui_->ions_list_widget->findItems("C-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> x_ions = ui_->ions_list_widget->findItems("X-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> y_ions = ui_->ions_list_widget->findItems("Y-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> z_ions = ui_->ions_list_widget->findItems("Z-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> pc_ions = ui_->ions_list_widget->findItems("Precursor", Qt::MatchFixedString);
      QList<QListWidgetItem*> nl_ions = ui_->ions_list_widget->findItems("Neutral losses", Qt::MatchFixedString);
      QList<QListWidgetItem*> ic_ions = ui_->ions_list_widget->findItems("Isotope clusters", Qt::MatchFixedString);
      QList<QListWidgetItem*> ai_ions = ui_->ions_list_widget->findItems("Abundant immonium-ions", Qt::MatchFixedString);

      OPENMS_PRECONDITION(a_ions.size() == 1, "String 'A-ions' does not exist in identification dialog.");
      OPENMS_PRECONDITION(b_ions.size() == 1, "String 'B-ions' does not exist in identification dialog.");
      OPENMS_PRECONDITION(c_ions.size() == 1, "String 'C-ions' does not exist in identification dialog.");
      OPENMS_PRECONDITION(x_ions.size() == 1, "String 'X-ions' does not exist in identification dialog.");
      OPENMS_PRECONDITION(y_ions.size() == 1, "String 'Y-ions' does not exist in identification dialog.");
      OPENMS_PRECONDITION(z_ions.size() == 1, "String 'Z-ions' does not exist in identification dialog.");
      OPENMS_PRECONDITION(pc_ions.size() == 1, "String 'Precursor' does not exist in identification dialog.");
      OPENMS_PRECONDITION(nl_ions.size() == 1, "String 'Neutral losses' does not exist in identification dialog.");
      OPENMS_PRECONDITION(ic_ions.size() == 1, "String 'Isotope clusters' does not exist in identification dialog.");
      OPENMS_PRECONDITION(ai_ions.size() == 1, "String 'Abundant immonium-ions' does not exist in identification dialog.");

      a_ions[0]->setCheckState(param_.getValue("preferences:idview:show_a_ions").toBool() ? Qt::Checked : Qt::Unchecked);
      b_ions[0]->setCheckState(param_.getValue("preferences:idview:show_b_ions").toBool() ? Qt::Checked : Qt::Unchecked);
      c_ions[0]->setCheckState(param_.getValue("preferences:idview:show_c_ions").toBool() ? Qt::Checked : Qt::Unchecked);
      x_ions[0]->setCheckState(param_.getValue("preferences:idview:show_x_ions").toBool() ? Qt::Checked : Qt::Unchecked);
      y_ions[0]->setCheckState(param_.getValue("preferences:idview:show_y_ions").toBool() ? Qt::Checked : Qt::Unchecked);
      z_ions[0]->setCheckState(param_.getValue("preferences:idview:show_z_ions").toBool() ? Qt::Checked : Qt::Unchecked);
      pc_ions[0]->setCheckState(param_.getValue("preferences:idview:show_precursor").toBool() ? Qt::Checked : Qt::Unchecked);
      nl_ions[0]->setCheckState(param_.getValue("preferences:idview:add_losses").toBool() ? Qt::Checked : Qt::Unchecked);
      ic_ions[0]->setCheckState(param_.getValue("preferences:idview:add_isotopes").toBool() ? Qt::Checked : Qt::Unchecked);
      ai_ions[0]->setCheckState(param_.getValue("preferences:idview:add_abundant_immonium_ions").toBool() ? Qt::Checked : Qt::Unchecked);
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
      p.setValue("preferences:default_path", ui_->default_path->text());
      p.setValue("preferences:default_path_current", ui_->default_path_current->isChecked());

      p.setValue("preferences:use_cached_ms1", ui_->use_cached_ms1->isChecked());
      p.setValue("preferences:use_cached_ms2", ui_->use_cached_ms2->isChecked());

      p.setValue("preferences:tmp_file_path", ui_->temp_path->text());
      p.setValue("preferences:number_of_recent_files", ui_->recent_files->value());
      p.setValue("preferences:default_map_view", ui_->map_default->currentText());
      p.setValue("preferences:intensity_cutoff", ui_->map_cutoff->currentText());
      p.setValue("preferences:on_file_change", ui_->on_file_change->currentText());

      p.setValue("preferences:1d:peak_color", ui_->color_1D->getColor().name());
      p.setValue("preferences:1d:highlighted_peak_color", ui_->selected_1D->getColor().name());
      p.setValue("preferences:1d:icon_color", ui_->icon_1D->getColor().name());

      p.setValue("preferences:2d:dot:gradient", ui_->peak_2D->gradient().toString());
      p.setValue("preferences:2d:mapping_of_mz_to", ui_->mapping_2D->currentText());
      p.setValue("preferences:2d:dot:feature_icon", ui_->feature_icon_2D->currentText());
      p.setValue("preferences:2d:dot:feature_icon_size", ui_->feature_icon_size_2D->value());

      p.setValue("preferences:3d:dot:gradient", ui_->peak_3D->gradient().toString());
      p.setValue("preferences:3d:dot:shade_mode", ui_->shade_3D->currentIndex());
      p.setValue("preferences:3d:dot:line_width", ui_->line_width_3D->value());

      // id view
      p.setValue("preferences:idview:a_intensity", ui_->a_intensity->value(), "Default intensity of a-ions");
      p.setValue("preferences:idview:b_intensity", ui_->b_intensity->value(), "Default intensity of b-ions");
      p.setValue("preferences:idview:c_intensity", ui_->c_intensity->value(), "Default intensity of c-ions");
      p.setValue("preferences:idview:x_intensity", ui_->x_intensity->value(), "Default intensity of x-ions");
      p.setValue("preferences:idview:y_intensity", ui_->y_intensity->value(), "Default intensity of y-ions");
      p.setValue("preferences:idview:z_intensity", ui_->z_intensity->value(), "Default intensity of z-ions");
      p.setValue("preferences:idview:relative_loss_intensity", ui_->relative_loss_intensity->value(), "Relativ loss in percent");
      p.setValue("preferences:idview:tolerance", ui_->tolerance->value(), "Alignment tolerance");

      QList<QListWidgetItem*> a_ions = ui_->ions_list_widget->findItems("A-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> b_ions = ui_->ions_list_widget->findItems("B-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> c_ions = ui_->ions_list_widget->findItems("C-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> x_ions = ui_->ions_list_widget->findItems("X-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> y_ions = ui_->ions_list_widget->findItems("Y-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> z_ions = ui_->ions_list_widget->findItems("Z-ions", Qt::MatchFixedString);
      QList<QListWidgetItem*> pc_ions = ui_->ions_list_widget->findItems("Precursor", Qt::MatchFixedString);
      QList<QListWidgetItem*> nl_ions = ui_->ions_list_widget->findItems("Neutral losses", Qt::MatchFixedString);
      QList<QListWidgetItem*> ic_ions = ui_->ions_list_widget->findItems("Isotope clusters", Qt::MatchFixedString);
      QList<QListWidgetItem*> ai_ions = ui_->ions_list_widget->findItems("Abundant immonium-ions", Qt::MatchFixedString);
      
      p.setValue("preferences:idview:show_a_ions", fromCheckState(a_ions[0]->checkState()), "Show a-ions");
      p.setValue("preferences:idview:show_b_ions", fromCheckState(b_ions[0]->checkState()), "Show b-ions");
      p.setValue("preferences:idview:show_c_ions", fromCheckState(c_ions[0]->checkState()), "Show c-ions");
      p.setValue("preferences:idview:show_x_ions", fromCheckState(x_ions[0]->checkState()), "Show x-ions");
      p.setValue("preferences:idview:show_y_ions", fromCheckState(y_ions[0]->checkState()), "Show y-ions");
      p.setValue("preferences:idview:show_z_ions", fromCheckState(z_ions[0]->checkState()), "Show z-ions");
      p.setValue("preferences:idview:show_precursor", fromCheckState(pc_ions[0]->checkState()), "Show precursor");
      p.setValue("preferences:idview:add_losses", fromCheckState(nl_ions[0]->checkState()), "Show neutral losses");
      p.setValue("preferences:idview:add_isotopes", fromCheckState(ic_ions[0]->checkState()), "Show isotopes");
      p.setValue("preferences:idview:add_abundant_immonium_ions", fromCheckState(ai_ions[0]->checkState()), "Show abundant immonium ions");

      if (param_.empty()) param_ = p; // if setParam() was not called before, param_ is empty and p is the only thing we have...
      else param_.update(p, false); // update with new values from 'p' to avoid loosing additional parameters and the existing descriptions already present in param_

      return param_;
    }

    void TOPPViewPrefDialog::browseDefaultPath_()
    {
      QString path = QFileDialog::getExistingDirectory(this, "Choose a directory", ui_->default_path->text());
      if (path != "")
      {
        ui_->default_path->setText(path);
      }
    }

    void TOPPViewPrefDialog::browseTempPath_()
    {
      QString path = QFileDialog::getExistingDirectory(this, "Choose a directory", ui_->temp_path->text());
      if (path != "")
      {
        ui_->temp_path->setText(path);
      }
    }

  }   //namespace Internal
} //namespace OpenMS

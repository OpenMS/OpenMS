/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you wish to add, delete or rename functions or slots use
** Qt Designer which will update this file, preserving your code. Create an
** init() function in place of a constructor, and a destroy() function in
** place of a destructor.
*****************************************************************************/

#include <qmessagebox.h>
#include <qmap.h>
#include <qstring.h>
#include <qsettings.h>
#include <qdir.h>
#include "SpecAnnotate.h"

void SettingsDialog::init()
{
    //rtti: does pa point to instance of SpecAnnotate?
    QWidget* pa = parentWidget();
    if (SpecAnnotate* pa_msa = dynamic_cast<SpecAnnotate*>(pa))
    {
	QMap<QString,QString>* settings;
	settings = pa_msa->getSettings();
	
	lineEdit1->setText((*settings)["db_username"]);
	lineEdit2->setText((*settings)["db_password"]);
	lineEdit3->setText((*settings)["db_host"]);
	lineEdit5->setText((*settings)["spl_path"]);	
	lineEdit6->setText((*settings)["peakfiles_path"]);
	lineEdit7->setText((*settings)["output_path"]);
    }	
}


void SettingsDialog::help()
{
QMessageBox::information(this, tr("Help: Settings Dialog"), 
			 tr("Please insert Settings valid for your system!"), 1);

}


void SettingsDialog::save()
{
    QSettings q_settings(QSettings::Ini);
    QDir current_dir;
    q_settings.insertSearchPath(QSettings::Unix, current_dir.absPath());
    
    q_settings.beginGroup("/specannotate/database"); 
    q_settings.writeEntry("db_username", lineEdit1->text());
    q_settings.writeEntry("db_password", lineEdit2->text());
    q_settings.writeEntry("db_host", lineEdit3->text());
    q_settings.endGroup();
    
    q_settings.beginGroup("/specannotate/paths");
    q_settings.writeEntry("spl_path", lineEdit5->text());
    q_settings.writeEntry("peakfiles_path", lineEdit6->text());
    q_settings.writeEntry("output_path", lineEdit7->text());

    q_settings.endGroup();
 
    //actualize settings in parent widget
    actualizeParentSettings();
    
}


void SettingsDialog::ok()
{
    save();
    accept();
}

void SettingsDialog::actualizeParentSettings()
{
    //rtti: does pa point to instance of SpecAnnotate?
    QWidget* pa = parentWidget();
    if (SpecAnnotate* pa_msa = dynamic_cast<SpecAnnotate*>(pa))
    {
	QMap<QString,QString>* settings;
	settings = pa_msa->getSettings();
	
	(*settings)["db_username"] = lineEdit1->text();
	(*settings)["db_password"] = lineEdit2->text();
	(*settings)["db_host"] = lineEdit3->text();
	(*settings)["spl_path"] = lineEdit5->text();
	(*settings)["peakfiles_path"] = lineEdit6->text();
	(*settings)["output_path"] = lineEdit7->text();
	
    }	
}


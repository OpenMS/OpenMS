/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you wish to add, delete or rename functions or slots use
** Qt Designer which will update this file, preserving your code. Create an
** init() function in place of a constructor, and a destroy() function in
** place of a destructor.
*****************************************************************************/

#include <qsqldatabase.h>
#include <qmap.h>
#include <qstring.h>
#include <qsqlquery.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include "SpecAnnotate.h"



void SampleDialog::init()
{
    //get settings from parent widget (rtti: does pa point to instance of SpecAnnotate?)
    QWidget* pa = this->parentWidget();
    if ((pa_msa = dynamic_cast<SpecAnnotate*>(pa)))
    {
	settings_ = pa_msa->getSettings();
    }
    else
    {
	exit(1);
    }

		//! the default database connection, using the QTDesigner-MySQL-Driver defined in OpenMS/include/OpenMS/config.h
		QSqlDatabase *defaultDB = QSqlDatabase::addDatabase( DB_PLUGIN );
    if ( ! defaultDB )
    {
	qWarning( "Failed to connect to driver" );
	
	((QStatusBar*)(pa_msa->statusBar()))->message( tr("Could not connect to Database"), 2000 );
    }
    defaultDB->setDatabaseName( "msannotate" );
    defaultDB->setUserName( (*settings_)["db_username"] );
    defaultDB->setPassword( (*settings_)["db_password"] );
    defaultDB->setHostName( (*settings_)["db_host"] );
    if ( ! defaultDB->open() )
    {
	qWarning( "Failed to open database: msannotate!" +
		  defaultDB->lastError().driverText() );
	qWarning( defaultDB->lastError().databaseText() );
	
	pa_msa->statusBar()->message( tr("Could not connect to Database"), 2000 );
    }
    
    //defaultDB is used automatically!!!!
    QSqlQuery query1( "SELECT identifier FROM protein;" );
    while ( query1.next() )
    {
	comboBox1->insertItem( query1.value( 0 ).toString());
    }	
    
    QSqlQuery query2( "SELECT enzyme_name FROM enzyme;" );
    while ( query2.next() )
    {
	comboBox2->insertItem( query2.value( 0 ).toString());
    }	
 
    
    QSqlQuery query3( "SELECT modification_name FROM modification ORDER BY modification_ID;" );
    while ( query3.next() )
    {
	listBox2->insertItem( query3.value( 0 ).toString());
    }	
     
    //insert possibility to select no enzyme
    comboBox2->insertItem("");
    
    //insert possible annotation methods
    comboBox3->insertItem("enumerate");
    comboBox3->insertItem("improved_enumerate");
    comboBox3->insertItem("peakwise_cormen");
    
    //insert known masstypes
    comboBox4->insertItem("average");
    comboBox4->insertItem("mono");
    
    //insert known peakfileformats
    comboBox5->insertItem("toll");
    comboBox5->insertItem("kerber");
    
    //load the default samplefile
    loadSample(((*settings_)["spl_path"] + "default.spl"));
}


QString SampleDialog::browse(QFileDialog::Mode mode, QString filetype)
{
    QString fn;
    if (mode == QFileDialog::ExistingFile)
    {
	if (filetype == "ini")
	{
	    fn = QFileDialog::getOpenFileName( (*settings_)["spl_path"], QString::null, this);
	}
	else if (filetype == "peak")
	{
	    fn = QFileDialog::getOpenFileName( (*settings_)["peakfiles_path"], QString::null, this);
	}
    }
    else if (mode == QFileDialog::DirectoryOnly)
    {
	fn = QFileDialog::getExistingDirectory( (*settings_)["output_path"], this);
    }
    else
    {
	fn = QString::null;
    }
    
    if ( !fn.isEmpty() )
    {
	return fn;
    }
    else
    {
	return QString::null;
    }
}	


void SampleDialog::browsePeakfile()
{
    QString fn = browse(QFileDialog::ExistingFile, "peak");
    if ((!fn.isEmpty()) && (!fn.isNull()))
    {
	lineEdit2->clear();
	lineEdit2->insert(fn);   
    }
}


void SampleDialog::browseOutputdir()
{
    QString fn = browse(QFileDialog::DirectoryOnly, "");
    if (	(!fn.isEmpty()) && (!fn.isNull()))
    {
	lineEdit3->clear();
	lineEdit3->insert(fn);   
    }
}


void SampleDialog::quit()
{
    pa_msa->close();
}


void SampleDialog::loadSample(QString filename)
{
    QString fn;
    if (filename == QString::null)
    {
	fn = browse(QFileDialog::ExistingFile, "ini");
	if ((!fn.isEmpty()) && (!fn.isNull()))
	{
	    lineEdit6->clear();
	    lineEdit6->insert(fn);   
	}
    }
    else
    {
	lineEdit6->clear();
	lineEdit6->insert(filename);   
	fn = filename;
    }
    
    //handle ini.file "by hand" with qt means
    QFile file( fn );
    if ( file.open( IO_ReadOnly ) ) 
    {
	QTextStream stream( &file );
	QString line;
	while ( !stream.atEnd() ) 
	{
	    line = stream.readLine(); 
	    if (line == "[SampleContents]")
	    {
		for (int i=0; i<2; i++)
		{
		    line = stream.readLine();
		    if (line.contains("enzyme="))
		    {
			line.remove("enzyme=");
			for (int i = 0; i < comboBox2->count(); i++)
			{	
			    if (comboBox2->text(i) == line)
			    {
				comboBox2->setCurrentItem(i);				    
			    }
			}
		    }
		    else if (line.contains("protein="))
		    {
			line.remove("protein=");
			for (int i = 0; i < comboBox1->count(); i++)
			{	
			    if (comboBox1->text(i) == line)
			    {
				comboBox1->setCurrentItem(i);				    
			    }
			}	
		    }
		}
	    }
	    else if (line == "[InputOutput]")
	    {
		for (int i=0; i<4; i++)
		{
		    line = stream.readLine();
		    if (line.contains("peakfile="))
		    {
			line.remove("peakfile=");
			lineEdit2->clear();
			lineEdit2->insert(line);
		    }
		    else if (line.contains("outputdir="))
		    {
			line.remove("outputdir=");
			lineEdit3->clear();
			lineEdit3->insert(line);
		    }
		    else if (line.contains("using_peakFile=true"))
		    {
			radioButton2->setChecked(true);
		    }
		    else if (line.contains("using_peakFile=false"))
		    {
			radioButton1->setChecked(true);
		    }
		    else if (line.contains("using_outputDir=true"))
		    {
			radioButton4->setChecked(true);
		    }
		    else if (line.contains("using_outputDir=false"))
		    {
			radioButton3->setChecked(true);
		    }
		}  
	    }
	    else if (line == "[Parameters]")
	    {
		for (int i=0; i<4; i++)
		{
		    line = stream.readLine();
		    if (line.contains("search_range="))
		    {
			line.remove("search_range=");
			lineEdit3_2->clear();
			lineEdit3_2->insert(line);
		    }
		    else if (line.contains("peakfile_format="))
		    {
			line.remove("peakfile_format=");
			for (int i = 0; i < comboBox5->count(); i++)
			{	
			    if (comboBox5->text(i) == line)
			    {
				comboBox5->setCurrentItem(i);				    
			    }
			}
		    }
		    else if (line.contains("masstype="))
		    {
			line.remove("masstype=");
			for (int i = 0; i < comboBox4->count(); i++)
			{	
			    if (comboBox4->text(i) == line)
			    {
				comboBox4->setCurrentItem(i);				    
			    }
			}
		    }
		    else if (line.contains("annotation_method="))
		    {
			line.remove("annotation_method=");
			for (int i = 0; i < comboBox3->count(); i++)
			{	
			    if (comboBox3->text(i) == line)
			    {
				comboBox3->setCurrentItem(i);				    
			    }
			}
		    }
		}  
	    }
	    else if (line == "[PartialModifications]")
	    {
		line = stream.readLine();
		textEdit1->clear();
		textEdit1->insert(line);
	    }
	    else if (line == "[OverallModifications]")
	    {
		listBox2->clearSelection();
		while (!stream.atEnd())
		{
		    line = stream.readLine();
		    for (uint i = 0; i < listBox2->count(); i++)
		    {	
			if (listBox2->text(i) == line)
			{
			    listBox2->setSelected(i, true);				    
			}
		    }
		}
	    }
	}
	file.close();
    }
}


void SampleDialog::saveSample()
{
    QFile file(lineEdit6->text());
    if ( file.open( IO_WriteOnly ) ) 
    {
	QTextStream stream( &file );
	stream << "[SampleContents]" << endl;
	stream << ("protein=" + comboBox1->currentText()) << endl;
	stream << ("enzyme=" + comboBox2->currentText()) << endl;
	
	stream << endl << endl;
	
	stream << "[InputOutput]" << endl;
	stream << ("peakfile=" + lineEdit2->text()) << endl;
	stream << ("outputdir=" + lineEdit3->text()) << endl;
	if (radioButton1->isChecked())
	{
	    stream << "using_peakFile=false" << endl;
	}
	else
	{
	    stream << "using_peakFile=true" << endl;
	}
	if (radioButton3->isChecked())
	{
	    stream << "using_outputDir=false" << endl;
	}
	else
	{
	    stream << "using_outputDir=true" << endl;
	}
	
	stream << endl << endl;
	
	stream << "[Parameters]" << endl;
	stream << ("search_range=" + lineEdit3_2->text()) << endl;
	stream << ("peakfile_format=" + comboBox5->currentText()) << endl;
	stream << ("masstype=" + comboBox4->currentText()) << endl;
	stream << ("annotation_method=" + comboBox3->currentText()) << endl;
		
	stream << endl << endl;
	
	stream << "[PartialModifications]" << endl;
	stream << textEdit1->text() << endl;
	
	stream << endl << endl;
	
	stream << "[OverallModifications]";
	for (uint i = 0; i < listBox2->count(); i++)
	{
	    if (listBox2->isSelected(i))
	    {
		stream << endl;
		stream << listBox2->text(i);
	    }
	}
        file.close();
    }
    
    pa_msa->statusBar()->message( tr("Sample " + lineEdit6->text() + " Saved!"), 2000 );
}


void SampleDialog::saveAs()
{
    QString fn = QFileDialog::getSaveFileName(QString::null, QString::null, this);
    if ((!fn.isEmpty()) && (!fn.isNull()))
    {
	lineEdit6->clear();
	lineEdit6->insert(fn);   
    }
    
    saveSample();
}


void SampleDialog::annotate()
{
    //check if wrong combination of input/output options is chosen
    if ((radioButton2->isChecked()) && (radioButton3->isChecked()))
    {
	QMessageBox::information(this, "Wrong selection", "Reading peaks from file and storing annotations as metadata in spectrum cannot be selected together. \nPlease correct your selection!");
    }
    else
    {
	__gnu_cxx::hash_map<std::string, std::string> sample_data;
	std::vector<std::string> ov_mods;
	sample_data["protein"] =  toStlString(comboBox1->currentText());
	sample_data["enzyme"] =  toStlString(comboBox2->currentText());
	sample_data["peakfile"] = toStlString(lineEdit2->text());
	
	//if annotations should be exported as metadata: take "" as outputdir
	if (radioButton3->isChecked())
	{
	    sample_data["outputdir"] = "";
	}
	else
	{
	    sample_data["outputdir"] = toStlString(lineEdit3->text());
	}	
	sample_data["search_range"] = toStlString(lineEdit3_2->text());
	sample_data["peakfile_format"] = toStlString(comboBox5->currentText());
	sample_data["masstype"] = toStlString(comboBox4->currentText());
	sample_data["annotation_method"] = toStlString(comboBox3->currentText());
	sample_data["partial_modification_string"] = toStlString(textEdit1->text());
	
	for (uint i = 0; i < listBox2->count(); i++)
	{	
	    if (listBox2->isSelected(i))
	    {	
		ov_mods.push_back(listBox2->text(i));
	    }
	}	
   
	//clear peaklist from peaks of previous runs
	peaklist.clear();
   
	//if corresponding Button is checked: get peaklist from spectrum
	if (radioButton1->isChecked())
	{
	    peaklist  = SpectrumMDIWindowEnhanced::getInstance()->getActiveSpectrumSelectedPeaks();
	    if (peaklist.empty())
	    {
		QMessageBox::information(this, "Missing Peaks", "Reading peaks from spectrum not successful, no peaks selected. \nPlease select peaks first!");
		return;
	    }
	}	
    
	Annotate* annotate = new Annotate(this, tr("Annotating Sample " + lineEdit6->text() + "..."));
	annotate->show();
	annotate->run(sample_data, peaklist, ov_mods, settings_);
    }
}


void SampleDialog::inputModifications()
{
 
   InputModifications* inputmod = new InputModifications(this, tr("Partial Modification Inpu"));
   inputmod->show();
}


QString SampleDialog::getProtein()
{
    return comboBox1->currentText();
}


int SampleDialog::getProteinSize()
{
    QSqlQuery query( "SELECT no_of_aminoacids FROM protein WHERE identifier = \"" +  comboBox1->currentText() + "\";" );
    int result = 0;
    if ( query.next() )
    {
	result = query.value( 0 ).toInt();
    }	
    return result;
}


void SampleDialog::insertPartialMod(QString mod)
{
    textEdit1->clear();
    textEdit1->insert(mod);
}


std::string SampleDialog::toStlString( QString s )
{
    std::ostringstream ofst;
    ofst << s;
    return ofst.str();
}


void SampleDialog::loadSampleNoDefault()
{
    loadSample(QString::null);
}


void SampleDialog::importPeaklistFromFile()
{
    lineEdit2->setEnabled(true);
    pushButton5->setEnabled(true);
    textLabel2_3->setEnabled(true);
    comboBox5->setEnabled(true);
    textLabel1_4->setEnabled(true);
    radioButton3->setEnabled(false);
    if (radioButton3->isChecked())
    {
	radioButton4->toggle();
    }
}


void SampleDialog::importPeaks()
{
    lineEdit2->setEnabled(false);
    textLabel2_3->setEnabled(false);
    comboBox5->setEnabled(false);
    textLabel1_4->setEnabled(false);
    pushButton5->setEnabled(false);
    radioButton3->setEnabled(true);
}


void SampleDialog::exportFiles()
{
    lineEdit3->setEnabled(true);
    textLabel2_2->setEnabled(true);
    pushButton5_2->setEnabled(true);
}


void SampleDialog::exportMetadata()
{
    lineEdit3->setEnabled(false);
    textLabel2_2->setEnabled(false);
    pushButton5_2->setEnabled(false);
}

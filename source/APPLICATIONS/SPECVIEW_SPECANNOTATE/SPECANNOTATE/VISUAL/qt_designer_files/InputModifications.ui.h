/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you wish to add, delete or rename functions or slots use
** Qt Designer which will update this file, preserving your code. Create an
** init() function in place of a constructor, and a destroy() function in
** place of a destructor.
*****************************************************************************/



void InputModifications::init()
{
    QWidget* pa = this->parentWidget();
    if ((sd = dynamic_cast<SampleDialog*>(pa)))
    {
	//get settings from parent widget (rtti: does pa point to instance of SpecAnnotate?)
	QWidget* pa_pa = sd->parentWidget();
	if ((msa = dynamic_cast<SpecAnnotate*>(pa_pa)))
	{
	    settings_ = msa->getSettings();
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
	    ((QStatusBar*)(msa->statusBar()))->message( tr("Could not connect to Database"), 2000 );
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
	    msa->statusBar()->message( tr("Could not connect to Database"), 2000 );
	}
	
	//positions
	QSqlQuery q1( "SELECT protein_ID FROM protein WHERE identifier = \"" + sd->getProtein() + "\";" );		QString prot_ID;	
	if ( q1.next() )
	{
	   prot_ID = q1.value( 0 ).toString();
	}
	
	listBox1->clear();
	textLabel1->setText("Positions in Protein " + sd->getProtein());
	QString insert;
	for (int i = 0; i < sd->getProteinSize(); i++)
	{
	    insert.setNum(i);	
	    QString aa_ID;
	    QSqlQuery q2( "SELECT aminoacid_ID FROM sequence WHERE protein_ID = " + prot_ID + 
			  " AND s_position = " + insert + ";" );   
	    if ( q2.next() )
	    {	
		aa_ID = q2.value( 0 ).toString();
	    }   
	    QSqlQuery q3( "SELECT three_letter_code FROM aminoacid WHERE aminoacid_ID = " + aa_ID + ";" );
	    if ( q3.next() )
	    {	
		insert += " (";
		insert += q3.value( 0 ).toString();
		insert += ")";
	    }   

	    
	    listBox1->insertItem(insert);
	 }	
	
	//modifications	
	listBox1_2->clear();
	QSqlQuery query( "SELECT modification_name FROM modification ORDER BY modification_ID;" );		while ( query.next() )
	{
	    listBox1_2->insertItem( query.value( 0 ).toString());
	}
	
    }
    else
    {
	exit(1);
    }
}



void InputModifications::done()
{
    //FUEGE TEXT DES VIEWS INS ENSPR FELD IN QSAMPLEDIALOG EIN
    QString mod = textBrowser1->text();
    mod.truncate(mod.length()-2);
    mod += "*";
    sd->insertPartialMod(mod);
    
    close();
}





void InputModifications::addGroup()
{
    QString group;
    std::list<int> int_group;
  
    //get all modifications of current selection (group)
    for (uint i = 0; i < listBox1_2->count(); i++)
	{	
	    if (listBox1_2->isSelected(i))
	    {
		QSqlQuery query( "SELECT modification_ID FROM modification WHERE modification_name = \"" + listBox1_2->text(i) + "\";" );
		if ( query.next() )
		{
		    int_group.push_back(query.value( 0 ).toInt());
		}
	    }
	}
    //to be sure same modification groups are recognized as same: sort
    int_group.sort();
    group = "( ";
    QString mod;
    for (std::list<int>::iterator it = int_group.begin(); it != int_group.end(); it++)
    {
	if (! (it == int_group.begin()))
	{
	    group += " , ";
	}
	mod.setNum(*it);
	group += mod;
    }
    group += " )";
 
    //get all positions of current selection (group)
    QString add_string; //this string will be actually added to the overall modification string
    bool is_first = true;
    for (uint i = 0; i < listBox1->count(); i++)
    {		
	if (listBox1->isSelected(i))
	{
	    {
		if (!is_first)
		{
		    add_string += " ; ";
		}
		add_string += (listBox1->text(i));
		//remove three-letter-codes
		add_string.remove((add_string.length() - 6), 6);
		add_string += " ";
		add_string += group; 
		
		is_first = false;
	    }
	}
    }
 
    add_string += " ; ";
    
    textBrowser1->insert(add_string);
    
    //clear selection
    resetSelection();
}


void InputModifications::resetSelection()
{
    listBox1->clearSelection();
    listBox1_2->clearSelection();
}


void InputModifications::resetString()
{
    textBrowser1->clear();
}

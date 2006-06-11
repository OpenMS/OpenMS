/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you wish to add, delete or rename functions or slots use
** Qt Designer which will update this file, preserving your code. Create an
** init() function in place of a constructor, and a destroy() function in
** place of a destructor.
*****************************************************************************/



void Annotate::run(__gnu_cxx::hash_map<std::string, std::string> sample_data, std::vector<OpenMS::Spectrum1DWidget::Spectrum1D::const_iterator>& peaklist, std::vector<std::string> ov_mods, QMap<QString, QString> * settings )
{
    //fill member variable
    settings_ = settings;
    setCaption("Annotating...");
     
    // for measuring the time
    t.start();                               //for counting time since start
    timer_ID = startTimer(1000);  //starting first one-second-timer for showing elapsed time every second
    db_display_update_timer = startTimer(60000); //first one-minute timer for updating db info	
    updateDBDisplay();
 
    //creating  a new Annotation-Thread (object is not deleted until window is closed)
    qathread = new AnnotateThread(sample_data, peaklist, ov_mods, (*settings)["db_username"] , (*settings)["db_password"] , (*settings)["db_host"], this);
    
    qathread->start();
}


void Annotate::addOutput(std::string s)
{
    textBrowser1->append(QString(s));
    update();
}


void Annotate::ready()
{
    killTimer(timer_ID);
    killTimer(db_display_update_timer);
    killTimers();
    textBrowser1->setContentsPos(0, 0);
    updateDBDisplay();
}


void Annotate::abort()
{
    //terminate annotation-thread
    if (qathread->running())
    {
	qathread->terminate();
	qathread->wait();
	killTimer(timer_ID);
	killTimer(db_display_update_timer);
	textBrowser1->setContentsPos(0, 0);
	updateDBDisplay();
	QMessageBox::information(this, tr("Warning:"), tr("Annotation aborted by user!"));
    }
}


void Annotate::closeWindow()
{
   //clean up
    if (qathread->running())
    {
	abort();
    }
    close();
}


void Annotate::timerEvent( QTimerEvent * e )
{
    if (e->timerId() == timer_ID)
    {
	//update lcd display
	QTime elapsed(0,0,0,0);
	elapsed = elapsed.addMSecs(t.elapsed());
	lCDNumber1->display(elapsed.toString()); //elapsed.toString("hh:mm:ss"));
	update();
	
	timer_ID = startTimer(1000);
	
    }
    else if (e->timerId() == db_display_update_timer)
    {
	updateDBDisplay();
    }
}


void Annotate::updateDBDisplay()
{
    //connect to database
    QSqlDatabase *defaultDB = QSqlDatabase::addDatabase( DB_PLUGIN );
    if ( ! defaultDB )
    {
	qWarning( "Failed to connect to driver" );
	return;	
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
	return;
    }
  
    //execute queries
    QSqlQuery annotations = defaultDB->exec("SELECT count(*) FROM annotation");
    QSqlQuery mod_comb = defaultDB->exec("SELECT count(*) FROM modification_combination");
    QSqlQuery real_mod = defaultDB->exec("SELECT count(*) FROM realized_modification");
    QSqlQuery mod_comb_posless = defaultDB->exec("SELECT count(*) FROM modification_combination_positionless");
    QSqlQuery real_mod_posless = defaultDB->exec("SELECT count(*) FROM realized_modification_positionless");
	
    //process them and update lcd dislplays 
    if (annotations.isActive()) 
    {
	annotations.next(); //first state is BEFORE value that we want
	lCDNumber2->display(annotations.value(0).toString());	
    }
    if (mod_comb.isActive()) 
    {	
	mod_comb.next(); //first state is BEFORE value that we want
	lCDNumber3->display(mod_comb.value(0).toString());	
    }
    if (real_mod.isActive()) 
    {
	real_mod.next(); //first state is BEFORE value that we want
	lCDNumber4->display(real_mod.value(0).toString());	
    }
    if (mod_comb_posless.isActive()) 
    {	
	mod_comb_posless.next(); //first state is BEFORE value that we want
	lCDNumber5->display(mod_comb_posless.value(0).toString());	
    }
    if (real_mod_posless.isActive()) 
    {
	real_mod_posless.next(); //first state is BEFORE value that we want
	lCDNumber6->display(real_mod_posless.value(0).toString());	
    }
    
    //schedule repaint
    update();
    
    //start new one-minute-timer for updating database info
    db_display_update_timer = startTimer(60000);
}


void Annotate::customEvent( QCustomEvent * e )
{
    if (e->type() == 65432)   //an OutputEvent
    {  
	OutputEvent* ue = (OutputEvent*) e;
	addOutput(ue->output());
    }
    else if (e->type() == 65433)  //a FinishEvent
    {
	ready();
	QMessageBox::information(this, tr("Notification:"), tr("Annotation of Peaks finished!"));
    }
}

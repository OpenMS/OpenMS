//! [Logger] 

  ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);  // set the log type (command line or a file)

    // set start progress (0) and end (ms_run.size() = the number of spectra)
    progresslogger.startProgress(0, ms_run.size(), "Doing some calculation...");

    for (PeakMap::Iterator it = ms_run.begin(); it != ms_run.end(); ++it)
    {
	  // update progress
         progresslogger.setProgress(ms_run.end() - it);
      
	  // do the actual calculations and processing ...
    }
    progresslogger.endProgress();

//! [Logger]

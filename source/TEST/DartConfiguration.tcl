# This file is configured by CMake automatically as DartConfiguration.tcl
# If you choose not to use CMake, this file may be hand configured, by
# filling in the required variables.


# Configuration directories and files
SourceDirectory: /Users/andreasbertsch/Data/OpenMS/OpenMS/source/TEST
BuildDirectory: /Users/andreasbertsch/Data/OpenMS/OpenMS/source/TEST

# Site is something like machine.domain, i.e. pragmatic.crd
Site: aldehyde.local

# Build name is osname-revision-compiler, i.e. Linux-2.4.2-2smp-c++
BuildName: Darwin-g++-4.2

# Submission information
IsCDash: TRUE
DropSite: www-bs2.informatik.uni-tuebingen.de
DropLocation: /services/OpenMS_services/CDash/submit.php?project=OpenMS
DropSiteUser: 
DropSitePassword: 
DropSiteMode: 
DropMethod: http
TriggerSite: http://www-bs2.informatik.uni-tuebingen.de/cgi-bin/Submit-Random-TestingResults.cgi
ScpCommand: /usr/bin/scp

# Dashboard start time
NightlyStartTime: 23:00:00 UTC

# Commands for the build/test/submit cycle
ConfigureCommand: "/Applications/CMake 2.6-4.app/Contents/bin/cmake" "/Users/andreasbertsch/Data/OpenMS/OpenMS/source/TEST"
MakeCommand: /opt/local/bin/gmake -i

# CVS options
# Default is "-d -P -A"
CVSCommand: /usr/bin/cvs
CVSUpdateOptions: -d -A -P

# Subversion options
SVNCommand: /usr/bin/svn
SVNUpdateOptions: 

# Generic update command
UpdateCommand: /usr/bin/svn
UpdateOptions: 
UpdateType: svn

# Dynamic analisys and coverage
PurifyCommand: 
ValgrindCommand: 
ValgrindCommandOptions: 
MemoryCheckCommand: MEMORYCHECK_COMMAND-NOTFOUND
MemoryCheckCommandOptions: 
MemoryCheckSuppressionFile: 
CoverageCommand: /usr/bin/gcov

# Testing options
# TimeOut is the amount of time in seconds to wait for processes
# to complete during testing.  After TimeOut seconds, the
# process will be summaily terminated.
# Currently set to 25 minutes
TimeOut: 1500

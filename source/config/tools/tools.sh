
# NOTE:  If you get errors when you say ". tools.sh" at a bash prompt,
# the reason might be that this file has DOS line endings.
# => then apply dos2unix 
# Have fun ;-)


## Clemens' stuff

chop()
{
	sed s/:/\\n/g -
}

# Usage e.g.:  echo $PATH | chop_verbose
chop_verbose() {
    echo '#BEGIN'
    chop | sed 's/^/\t/' -
    echo '#END'
}

function path () {
    echo $PATH | chop
}

function ldli () {
    echo $LD_LIBRARY_PATH | chop
}

function g () {
    echo '###' gvim --remote $*
    gvim --remote $*
}

function i ()
{
    echo You are `whoami` on `hostname`.
}

# Show disk usage
function dusn ()
{
    du $* | sort -n
}

# Change to "real" working directory
function cwd ()
{
	cd "`pwd -P`i"	
}

# go to .../source/...
function gos()
{
	cd `pwd | sed "s/include\/OpenMS/source/"`
}

# go to .../include/OpenMS/...
function goi()
{
	cd `pwd | sed "s/source/include\/OpenMS/"`
}

# go to .../OpenMS/...
function godev()
{
	cd `pwd | sed "s:Release:OpenMS:"` ; ldli | grep --color Release
}

# go to .../Release/...
function gorel()
{
	cd `pwd | sed "s:OpenMS:Release:"` ; ldli | grep --color OpenMS
}


###  Special path to QT4  (used if present)
export OPENMS_QT4=${HOME}/Qt4

###  Location of OpenMS contrib package  (used if present)
export OPENMS_CONTRIB=${HOME}/contrib

function prefix_PATH ()
{
    export PATH=$1:$PATH
    #echo $PATH | chop_verbose
}

function prefix_PATH. () {
    export PATH=`pwd`:$PATH
    #echo $PATH | chop_verbose
}

function postfix_PATH () {
    export PATH=$PATH:$1
    #echo $PATH | chop_verbose
}

function prefix_LD_LIBRARY_PATH () {
    export LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
    #echo $LD_LIBRARY_PATH | chop_verbose
}

function prefix_LD_LIBRARY_PATH. () {
    export LD_LIBRARY_PATH=`pwd`:$LD_LIBRARY_PATH
    #echo $LD_LIBRARY_PATH | chop_verbose
}

function postfix_LD_LIBRARY_PATH () {
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$1
    #echo $LD_LIBRARY_PATH | chop_verbose
}

function OpenMS-DIR ()
{
    if [ ! -d $1"" ]; then
	echo \'$1\' " is not a directory"
    else
	echo "Setting up environment for OpenMS in directory " $1
	echo
	export OPENMS=$1
	echo export OPENMS=${OPENMS}
	export PROJECT=$1
	echo export PROJECT=${PROJECT}
	export CXX='g++'
	echo export CXX=\'${CXX}\'
	export CXXCPP='g++ -E'
	echo export CXXCPP=\'${CXXCPP}\'
	export CC=gcc
	echo export CC=\'${CC}\'
	export CPP='gcc -E'
	echo export CPP=\'${CPP}\'
	echo
	if [ -d ${OPENMS_CONTRIB}"" ]; then
	    echo "### Adding " ${OPENMS_CONTRIB} " to LD_LIBRARY_PATH"
    	    prefix_LD_LIBRARY_PATH ${OPENMS_CONTRIB}
	fi
	if [ -d ${OPENMS_QT4}"" ]; then
	    echo "### Adding " ${OPENMS_QT4}/lib " to LD_LIBRARY_PATH"
	    prefix_LD_LIBRARY_PATH ${OPENMS_QT4}/lib
	    echo "### Adding " ${OPENMS_QT4}/bin " to PATH"
	    prefix_PATH ${OPENMS_QT4}/bin
	fi
	echo "### Adding " ${OPENMS}/lib " to LD_LIBRARY_PATH"
    	prefix_LD_LIBRARY_PATH ${OPENMS}/lib
	echo "### Adding " ${OPENMS}/bin " to PATH"
    	prefix_PATH ${OPENMS}/bin
    	echo
	echo cd ${OPENMS}
	cd ${OPENMS}
    fi
}

function OpenMS-develop ()
{
	OpenMS-DIR ${HOME}/OpenMS
}
# I tend to misspell this important command, so here is a handy shortcut
alias om=OpenMS-develop

function OpenMS-release ()
{
	OpenMS-DIR ${HOME}/Release
}


export PATCH_VERSION_CONTROL=numbered
export VERSION_CONTROL=numbered
export GREP_OPTIONS='--color'

function pl () {
export PROMPT_LENGTH=$1
}
pl 20

function reset-prompt () {
echo reset-prompt
case $TERM in
	xterm*|rxvt*)
		TITLEBAR='\[\033]0;\u@\h:\w\007\]'
		;;
	*)
		TITLEBAR=""
		;;
esac;
case $TERM in
	dumb*)
		export PS1="\$(lPWD=\${#PWD}; len=\$PROMPT_LENGTH; len2=len-6;\
		if [ \${lPWD} -gt \${len} ]; then echo \"...\${PWD:\$((\${#PWD}-len2))} \";\
		else echo \\w\ ; fi)\$ "
		;;
 *)
		export PS1="${TITLEBAR}\[\033[0m\033[0;1m\]\$(lPWD=\${#PWD}; len=\$PROMPT_LENGTH; len2=len-6;\
 		if [ \${lPWD} -gt \${len} ]; then echo \"\[\033[0;31m\]...\[\033[m\033[0;1m\]\${PWD:\$((\${#PWD}-len2))} \";\
 		else echo \\w\ ; fi)\[\033[0;32m\]\$ \[\033[0;0m\]"
		;;
 esac;
}

function P () {
reset-prompt
export PS1=$PS1"\[\e]30;\H:\w\a\]"
}

function p () {
reset-prompt
}

case $TERM in
	dumb*|eterm*)
		p
		;;
	*)
		P
		;;
esac;



### Marc's stuff below

### 
### All the scripts assume that your OpenMS installation path
### is stored in the environment variable $OMS!
###

### readable svn diff output
alias svndiff="svn diff --diff-cmd diff -x -uwB | less"

### update OpenMS and show the status in the command line
u() { cd `echo $OMS`; svn update; rm -rf source/TEST/TOPP/Inspect_FASTAFile_test4.tmp.trie source/TEST/TOPP/Inspect_FASTAFile_test4.tmp.index; svn status | sort; }

### build and execute an OpenMS test
t() { cd `echo $OMS`/source/TEST; make $1_test && $1_test -V; }
### build and execute a TOPP test
tt() { cd `echo $OMS`/source/APPLICATIONS/TOPP/ && make $1 && cd ../../TEST/TOPP/ && make DEBUG=1 VERBOSE=1 $1_test; }


### list all C++ source files of OpenMS
sources() { find $(echo $OMS) -name "*.h" -o -name "*.C"; }
### find a string in all C++ source files of OpenMS
fsources() { sources | xargs grep $@; }
### replace a string in all C++ source files of OpenMS
rsources() { sources | xargs sed -i "s/$1/$2/g"; }

### builds and installs the lib
mi() { cd `echo $OMS`/source && make && make install; }
### builds, installs and tests the lib
mit() { mi && make test; }
### builds, installs and tests the lib and TOPP
mitt() { mit && make TOPP && make TOPPtest; }

# This is a bash.
echo This is a $0.
echo OpenMS environment set up.

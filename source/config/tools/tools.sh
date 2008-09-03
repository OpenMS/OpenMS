
# NOTE:  If you get errors when you say ". tools.sh" at a bash prompt,
# the reason might be that this file has DOS line endings.
# => then apply dos2unix 
# Have fun ;-)


## Clemens' stuff



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

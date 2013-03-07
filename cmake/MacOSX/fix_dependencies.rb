#!/usr/bin/ruby

# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

require 'getoptlong'
require 'pathname'
require 'set'
require 'fileutils'

##### GLOBAL VARIABLES
$handledLibraries = Set.new
lib = nil
bin = nil
$install_name_tool=`which install_name_tool`.strip
$DEBUG = false
$currentIndent=0
$executableId="@executable_path/../lib/"

###############################################################################
def debug(message)
  if($DEBUG)
    prefix=""
    for index in 0 ... $currentIndent
      prefix += "--"
    end
    puts "#{prefix}> #{message}"
  end
end

###############################################################################
def fixable(path)
  if path.match(/^\./)
    return false
  elsif path.match(/\.app$/)
    return false
  else
    return true
  end
end

###############################################################################
def isFramework(path)
  return path.match(/\.framework/)
end

###############################################################################
def cleanOtoolEntry(otool_line)
  return `echo "#{otool_line}" |  tr -d ' \t' | sed 's/(.*//g'`.strip
end

###############################################################################
def getId(otool_output)
  id = otool_output.delete_at(0).strip
  id = otool_output.delete_at(0).strip
  # clean id
  id = cleanOtoolEntry(id)
#  debug "Found ID: #{id}"  
  return id, otool_output
end

##### METHODS USED TO HANDLE THE BUNDLING

###############################################################################
def fixLoadPath(oldpath,libpath,target)
  debug "#{$install_name_tool} -change #{oldpath} #{$executableId}#{libpath} #{target}"
  `#{$install_name_tool} -change #{oldpath} #{$executableId}#{libpath} #{target}`
end

###############################################################################
def fixId(target, libname)
  debug "#{$install_name_tool} -id #{$executableId}#{libname} #{target}"
  `#{$install_name_tool} -id #{$executableId}#{libname} #{target}`
end

###############################################################################
def handleDependencies(otool_out, targetPath, currentLib)
  for index in 0 ... otool_out.size
    fix_lib=cleanOtoolEntry(otool_out[index])

    # (\/usr\/lib|\/System)
    if fix_lib.match(/^(\/usr\/lib|\/System)/)
      debug "Ignoring system-lib: #{fix_lib}"
    elsif fix_lib.match(/^@executable/)
      puts "Ignoring libs that were already relinked (#{fix_lib})"
    elsif not fix_lib.match(/\//) # we need a path here, otherwise it is a lib in the same directory
      debug "we do not fix libs which are in the same directory #{fix_lib}"
    else
      newPath=""
      libname=""
      if isFramework(fix_lib)
        newPath, libname = handleFramework(fix_lib, targetPath)
      else
        newPath, libname = handleDyLib(fix_lib, targetPath)
      end
      
      # fix loading of this library
      fixLoadPath(fix_lib, libname, currentLib)
    end
  end
end

###############################################################################
def getDylibName(lib)
  return File.basename(lib)
end

###############################################################################
def copyLib(lib, targetPath)
  libBasename=File.basename(lib)
  # new path
  newPath="#{targetPath}/#{libBasename}"
  
  if not File.exist?(newPath)
    debug "Copy #{libBasename} from #{lib} to #{targetPath}"
    `cp #{lib} #{targetPath}`
    `chmod u+rw #{newPath}`
  else
    debug "#{libBasename} exists in #{targetPath}"
  end
    
  return newPath, libBasename
end

###############################################################################
def getFrameworkName(frameworkPath)
  # get framework path (name of dir containing .framework)
  frameworkDir=frameworkPath.gsub(/(.*\.framework).*/,'\1')
  frameworkName=frameworkDir.gsub(/(.*\/)([^\/]*\.framework)/,'\2')
  return frameworkName  
end

###############################################################################
def getFrameworkDir(frameworkPath)
  frameworkDir=frameworkPath.gsub(/(.*\.framework).*/,'\1')
  return frameworkDir
end

###############################################################################
def copyFramework(frameworkPath, targetPath)
  # get framework path (name of dir containing .framework)
  frameworkDir=frameworkPath.gsub(/(.*\.framework).*/,'\1')
  frameworkName=frameworkDir.gsub(/(.*\/)([^\/]*\.framework)/,'\2')

  # the actual lib name
  internFrameworkDir=frameworkPath.gsub(/(.*\.framework)\/(.*)/,'\2')  
  libname=frameworkName+'/'+internFrameworkDir
  # the new path
  newFrameworkPath="#{targetPath}/#{libname}"

  if not File.exists?(targetPath + frameworkName)
    debug "Copy fw #{frameworkName} from #{frameworkDir} to #{targetPath}"
    `cp -r #{frameworkDir} #{targetPath}`
    # adjust rights
    `chmod u+rw #{newFrameworkPath}`
  else
    debug "fw #{frameworkName} already exists in #{targetPath}"
  end
    
  return newFrameworkPath, libname
end

###############################################################################
def handleFramework(frameworkPath, targetPath)
  $currentIndent+=1

  # copy framework to target directory
  newFrameWorkPath, libname=copyFramework(frameworkPath, targetPath)

  if not $handledLibraries.include?(libname)
    # run otool
    otool_out=`otool -L #{frameworkPath}`.strip.split(/\n/)
    id, otool_out = getId(otool_out)
    debug "Handle FW #{frameworkPath} -> #{id}"
  
    # update the id
    fixId(newFrameWorkPath,libname)

    # check the actual dependencies
    handleDependencies(otool_out, targetPath, newFrameWorkPath)
  
    # mark as processed
    $handledLibraries.add(libname)
  else
    debug "Already fixed #{libname}"
  end
  
  # update ids of dependencies
  $currentIndent-=1
  
  return newFrameWorkPath, libname
end

###############################################################################
def handleDyLib(dylibPath, targetPath)
  $currentIndent+=1

  # copy if necessary
  newDyLibPath,libname=copyLib(dylibPath, targetPath)
  
  if not $handledLibraries.include?(libname)
    # run otool
    otool_out=`otool -L #{dylibPath}`.strip.split(/\n/)
    id, otool_out = getId(otool_out)
    debug "Handle DYLIB #{dylibPath} -> #{id}"
  
    # fixId
    fixId(newDyLibPath,libname)
  
    # check the actual dependencies
    handleDependencies(otool_out, targetPath, newDyLibPath)
  
    # mark as processed
    $handledLibraries.add(libname)
  else
    debug "Already fixed #{libname}"
  end
  
  # readjust
  $currentIndent-=1
  
  return newDyLibPath,libname
end

###############################################################################
def handleBinary(binaryPath, targetPath)
  $currentIndent+=1
  
  debug "Fixing #{binaryPath}"
  
  # no copy, no id change; juts run otool
  otool_out=`otool -L #{binaryPath}`.strip.split(/\n/)

  # remove 1st line containing the binary name
  otool_out.delete_at(0) 

  # handle all referenced libraries and copy them if necessary to 
  # the lib ref path
  handleDependencies(otool_out, targetPath, binaryPath)
  
  # readjust
  $currentIndent-=1
end

###############################################################################
###############################################################################
###############################################################################

##### INSTANTIATE CMD PARSER
opts = GetoptLong.new(
  [ '--help', '-h', GetoptLong::NO_ARGUMENT ],
  [ '--verbose', '-v', GetoptLong::NO_ARGUMENT ],
  [ '--lib-path', '-l', GetoptLong::REQUIRED_ARGUMENT ],
  [ '--install-name-tool', '-i', GetoptLong::OPTIONAL_ARGUMENT ],
  [ '--bin-path', '-b', GetoptLong::REQUIRED_ARGUMENT ]
)

usage = "#{File.basename($0)} --bin-path PATH-TO-BINARIES --lib-path PATH-TO-STORE-LIBRARIES
-h, --help:
  show help
-i, --install-name-tool:
  specify the install-name-tool manually
-l, --lib-path:
  the target path where the handled libraries should be stored and possibly already libraries exist
-b, --bin-path: 
  the path were all the binaries that should be handled are located
-v, --verbose:
  increase verbosity
"

opts.each do |opt, arg|
  case opt
    when '--help'
      puts "#{usage}"
    when '--install-name-tool'
      $install_name_tool = arg
      puts "Use #{install_name_tool} to fix binaries in #{bin}"
    when '--lib-path'
      lib = Pathname.new(arg).realpath
    when '--bin-path'
      bin = Pathname.new(arg).realpath
    when '-v'
      $DEBUG = true
  end
end

if lib == nil or bin == nil
  puts "Please provide a bin and lib path"
  puts "#{usage}"
  exit 1
end

# fix libraries contained in lib-path
for content in Dir.entries(lib) 
  if fixable(content)
    if isFramework(content)
#      handleFramework(lib + content, lib)
    else
      handleDyLib(lib + content, lib)
    end
  end
end

# fix binary references
for content in Dir.entries(bin)
  if fixable(content)
    handleBinary(bin + content, lib)
  end
end
  

#!/usr/bin/ruby

# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
$lib_dir = nil
$bin_dir = nil
$plugin_dir = nil
$install_name_tool=`which install_name_tool`.strip
$DEBUG = false
$currentIndent=0
$executableId="@executable_path/"
$EXTRACTFW = false

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
def fixable(name, path)
  if File.directory?(path + name)
    return false
  elsif name.match(/^\./)
    return false
  elsif name.match(/\.app$/)
    return false
  else
    filename = "#{path}/#{name}"
    debug filename.to_s
    otool_out=`otool -L #{filename}`
    ## Otool returns exit code 0 even if it is not an object file -> check output
    return !( otool_out =~ /.*object file.*/ )
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
def extractInstallName(otool_output)
  # Skip file basename
  otool_output.delete_at(0).strip
  # Save install name (id)
  id = otool_output.delete_at(0).strip
  # clean id
  id = cleanOtoolEntry(id)
  # return id and remaining otool output = dependencies
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
    elsif fix_lib.start_with?("@")
      puts "Ignoring libs that are referenced from a relative reference (#{fix_lib})"
      if not fix_lib.start_with?($executableId)
        puts "Warning: (#{fix_lib}) does not match the requested prefix, though."
      end
    elsif not fix_lib.match(/\//) # we need a path here, otherwise it is a lib in the same directory
      debug "we do not fix libs which are in the same directory #{fix_lib}"
    else
      newPath=""
      libname=""
      if isFramework(fix_lib)
        _, libname = handleFramework(fix_lib, $lib_dir)
      else
        _, libname = handleDyLib(fix_lib, $lib_dir)
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
  frameworkDir=frameworkPath.to_s.gsub(/(.*\.framework).*/,'\1')
  frameworkName=frameworkDir.to_s.gsub(/(.*\/)([^\/]*\.framework)/,'\2')

  # the actual lib name
  internFrameworkDir=frameworkPath.to_s.gsub(/(.*\.framework)\/(.*)/,'\2')
  libname=frameworkName+'/'+internFrameworkDir
  # the new path
  newFrameworkPath="#{targetPath}/#{libname}"

  if not File.exist?(targetPath + frameworkName)
    debug "Copy fw #{frameworkName} from #{frameworkDir} to #{targetPath}"
    # preserve symlinks
    `cp -R #{frameworkDir} #{targetPath}`
    # adjust rights
    `chmod -R u+rwX #{newFrameworkPath}`
  else
    debug "fw #{frameworkName} already exists in #{targetPath}"
  end

  return newFrameworkPath, libname
end

###############################################################################
def copyLibFromFramework(frameworkPath, targetPath)
  libname=File.basename(frameworkPath)

  # the new path
  newFrameworkPath="#{targetPath}/#{libname}"

  if not File.exist?(targetPath + libname)
    debug "Copy fw #{frameworkName} from #{frameworkDir} to #{targetPath}"
    # preserve symlinks
    `cp #{frameworkPath} #{newFrameworkPath}`
    # adjust rights
    `chmod u+rwX #{newFrameworkPath}`
  else
    debug "fw #{frameworkName} already exists in #{targetPath}"
  end

  return newFrameworkPath, libname
end

###############################################################################
def handleFramework(frameworkPath, targetPath)
  $currentIndent+=1

  if $EXTRACTFW
    newFrameWorkPath, libname = copyLibFromFramework(frameworkPath, targetPath)
  else
    # copy framework to target directory
    newFrameWorkPath, libname= copyFramework(frameworkPath, targetPath)
  end

  if not $handledLibraries.include?(libname)
    # run otool
    otool_out=`otool -L #{frameworkPath}`.strip.split(/\n/)
    # strips first two lines
    id, otool_out = extractInstallName(otool_out)

    # update the install_name (-id)
    fixId(newFrameWorkPath,libname)
    debug "Handle FW #{frameworkPath} -> #{id}"

    # check the actual dependencies (lines 3++)
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
    # strips first two lines
    id, otool_out = extractInstallName(otool_out)

    # update install_name (-id) of current DyLib
    fixId(newDyLibPath,libname)
    debug "Handle DYLIB #{dylibPath} --> #{id}"

    # check the actual dependencies (lines 3++)
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
def handleBinary(binaryPath)
  $currentIndent+=1

  debug "Fixing binary #{binaryPath}"

  # no copy, no id change; juts run otool
  otool_out=`otool -L #{binaryPath}`.strip.split(/\n/)

  # remove 1st line containing the binary name
  otool_out.delete_at(0)

  # handle all referenced libraries and copy them if necessary to
  # the lib ref path
  handleDependencies(otool_out, $lib_dir, binaryPath)

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
  [ '--bin-path', '-b', GetoptLong::OPTIONAL_ARGUMENT ],
  [ '--plugin-path', '-p', GetoptLong::OPTIONAL_ARGUMENT ],
  [ '--path-prefix', '-e', GetoptLong::OPTIONAL_ARGUMENT ],
  [ '--extract-from-framework', '-f', GetoptLong::NO_ARGUMENT]
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
-p, --plugin-path:
  the path were optional plugin libraries like from QT5 are located
-e, --path-prefix:
  the prefix that is added to the new install_name (default: @executable_path/)
-f, --extract-from-framework:
  extract the linked libraries from their Framework folder. CAUTION: this does not copy Headers and Resources. Not tested with already present frameworks.
-v, --verbose:
  increase verbosity
"

opts.each do |opt, arg|
  case opt
    when '--help'
      puts usage.to_s
    when '--install-name-tool'
      $install_name_tool = arg
      puts "Use #{install_name_tool} to fix binaries in #{bin}"
    when '--lib-path'
      $lib_dir = Pathname.new(arg).realpath
    when '--bin-path'
      $bin_dir = Pathname.new(arg).realpath
    when '--plugin-path'
      $plugin_dir = Pathname.new(arg).realpath
    when '--path-prefix'
      $executableId = arg
    when '--extract-from-framework'
      $EXTRACTFW = true
    when '--verbose'
      $DEBUG = true
  end
end

if $lib_dir.nil?
  puts "Please provide at least a lib path"
  puts usage.to_s
  exit 1
elsif !$bin_dir.nil?
  $executableId = $executableId + $lib_dir.relative_path_from($bin_dir).to_s
  $executableId += "/"
  puts "Substituting prefix to find libs with:"
  puts $executableId
end

debug "HANDLING LIB DIR"
# fix libraries contained in lib-path
# recurse two-levels to capture the libraries inside the frameworks
Dir.chdir($lib_dir.to_s) do
  lib_files = Dir.glob(["*","*/*"])
  for content in lib_files
    if fixable(content, $lib_dir)
      if isFramework(content)
        handleFramework($lib_dir + content, $lib_dir)
      elsif (content.end_with?(".dylib") or content.end_with?(".so"))
        handleDyLib($lib_dir + content, $lib_dir)
      else
        debug "Skipped #{$lib_dir + content} -- No lib or framework?"
      end
    else
      debug "Skipped #{$lib_dir + content} -- Otool not executable on it."
    end
  end
end

debug "HANDLING BIN DIR"
# fix binary references
if !$bin_dir.nil?
  for content in Dir.entries($bin_dir)
    if fixable(content, $bin_dir)
      if (content.end_with?(".dylib") or content.end_with?(".so"))
        debug "Handle dylib #{$bin_dir + content}"
        handleDyLib($bin_dir + content, $bin_dir)
      else
        debug "Handle binary #{$bin_dir + content}"
        handleBinary($bin_dir + content)
      end
    else
      debug "Skipped #{$bin_dir + content} -- No binary or object?"
    end
  end
end

if !$plugin_dir.nil?
  debug "HANDLING PLUGIN DIR"
  for content in Dir.glob("#{$plugin_dir}/**/*.dylib")
    debug "Handle dylib #{$plugin_dir + content}"
    handleDyLib($plugin_dir + content, File.dirname($plugin_dir + content))
  end
end


#!/usr/bin/ruby

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
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
$auto_relative = true
$nocopy = false

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
def handleDependencies(otool_out, targetPath, currentLib, rpaths)
  rpath = %x[otool -l #{currentLib} | grep -A2 LC_RPATH | grep path | sed -n "s/^.*path\\s*\\(.*\\)(offset.*\$/\\1/p"]
  rpaths += rpath.split(/\n/)

  for index in 0 ... otool_out.size
    fix_lib=cleanOtoolEntry(otool_out[index])

    # (\/usr\/lib|\/System)
    if fix_lib.match(/^(\/usr\/lib|\/System)/)
      debug "Ignoring system-lib: #{fix_lib}"
    elsif $nocopy
      if isFramework(fix_lib)
        # get framework path (name of dir containing .framework)
        frameworkDir=fix_lib.to_s.gsub(/(.*\.framework).*/,'\1')
        frameworkName=frameworkDir.to_s.gsub(/(.*\/)([^\/]*\.framework)/,'\2')

        # the actual lib name
        internFrameworkDir=fix_lib.to_s.gsub(/(.*\.framework)\/(.*)/,'\2')
        libname=frameworkName+'/'+internFrameworkDir
        fixLoadPath(fix_lib, libname, currentLib)
      else
        fixLoadPath(fix_lib, File.basename(fix_lib), currentLib)
      end
    elsif fix_lib.start_with?("@loader") or fix_lib.start_with?("@executable")
      if $EXTRACTFW
        # only fix loading of this library if it does not start with executable/loader_path.
        # It basically means it should be fixed already. If it was not we will not be able to
        # find it anyway.
        fixLoadPath(fix_lib, File.basename(fix_lib), currentLib)
      else
        debug "Ignoring libs that are referenced from a relative reference (#{fix_lib})"
        if not fix_lib.start_with?($executableId)
          puts "Warning: (#{fix_lib}) does not match the requested prefix, though."
        end
      end
    elsif not fix_lib.match(/\//) # we need a path here, otherwise it is a lib in the same directory
      # TODO combine with previous if-case? this is basically relative
      debug "We do not fix libs which are in the same directory: #{fix_lib}"
    else
      newPath=""
      libname=""
      if isFramework(fix_lib)
        _, libname = handleFramework(fix_lib, $lib_dir, rpaths)
      else
        _, libname = handleDyLib(fix_lib, $lib_dir, rpaths)
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

  if not File.exist?(newPath) and not lib.to_s.start_with?($lib_dir.to_s)
    debug "Copy #{libBasename} from #{lib} to #{targetPath}"
    `cp #{lib} #{targetPath}`
    `chmod u+rw #{newPath}`
  else
    debug "#{libBasename} exists in #{targetPath} or subfolder"
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

  if not File.exist?(targetPath + frameworkName) and not frameworkPath.to_s.start_with?($lib_dir.to_s)
    debug "Copy fw #{frameworkName} from #{frameworkDir} to #{targetPath}"
    # preserve symlinks
    `cp -R #{frameworkDir} #{targetPath}`
    # adjust rights
    `chmod -R u+rwX #{newFrameworkPath}`
  else
    debug "fw #{frameworkName} already exists in #{targetPath} or subfolder"
  end

  return newFrameworkPath, libname
end

###############################################################################
def copyLibFromFramework(frameworkPath, targetPath)
  libname=File.basename(frameworkPath)

  # the new path
  newFrameworkPath="#{targetPath}/#{libname}"

  if not File.exist?(targetPath + libname) and not frameworkPath.to_s.start_with?($lib_dir.to_s)
    debug "Copy lib #{libname} from #{frameworkPath} to #{targetPath}"
    # preserve symlinks
    `cp #{frameworkPath} #{newFrameworkPath}`
    # adjust rights
    `chmod u+rwX #{newFrameworkPath}`
  else
    debug "lib #{libname} already exists in #{targetPath} or subfolder"
  end

  return newFrameworkPath, libname
end

###############################################################################
def handleFramework(frameworkPath, targetPath, rpaths)
  $currentIndent+=1

  debug "RPATHS: #{rpaths}"

  if frameworkPath.to_s.start_with?("@rpath")
    for index in 0 ... rpaths.size
      if File.file?(frameworkPath.gsub("@rpath", rpaths[index].strip))
        frameworkPath = frameworkPath.gsub("@rpath", rpaths[index].strip)
        debug "Found lib with rpath at #{frameworkPath}"
        break
      end
    end
  end

  if not frameworkPath.to_s.start_with?("@rpath")
      if $EXTRACTFW
        newFrameWorkPath, libname = copyLibFromFramework(frameworkPath, targetPath)
      else
        # copy framework to target directory
        newFrameWorkPath, libname = copyFramework(frameworkPath, targetPath)
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
        handleDependencies(otool_out, targetPath, newFrameWorkPath, rpaths)

        # mark as processed
        $handledLibraries.add(libname)
      else
        debug "Already fixed #{libname}"
      end
  else
    debug "RPath but no lib found. Assuming fixed already: #{frameworkPath}"
    if !$EXTRACTFW
      # But at least adjust the libname which is used in handleDependencies to fix loading.
      # get framework path (name of dir containing .framework)
      frameworkDir=frameworkPath.to_s.gsub(/(.*\.framework).*/,'\1')
      frameworkName=frameworkDir.to_s.gsub(/(.*\/)([^\/]*\.framework)/,'\2')

      # the actual lib name
      internFrameworkDir=frameworkPath.to_s.gsub(/(.*\.framework)\/(.*)/,'\2')
      libname=frameworkName+'/'+internFrameworkDir
    else
      libname=File.basename(frameworkPath)
    end
    # the new path
    newFrameworkPath="#{targetPath}/#{libname}"
  end

  # update ids of dependencies
  $currentIndent-=1

  return newFrameWorkPath, libname
end

###############################################################################
def handleDyLib(dylibPath, targetPath, rpaths)
  $currentIndent+=1

  debug "RPATHS: #{rpaths}"

  if dylibPath.to_s.start_with?("@rpath")
    for index in 0 ... rpaths.size
      if File.file?(dylibPath.gsub("@rpath", rpaths[index].strip))
        dylibPath = dylibPath.gsub("@rpath", rpaths[index].strip)
        debug "Found lib with rpath at #{dylibPath}"
        break
      end
    end
  end

  if not dylibPath.to_s.start_with?("@rpath")
   # copy if necessary
   newDyLibPath,libname=copyLib(dylibPath, targetPath)

   if not $handledLibraries.include?(libname)
     # run otool
     otoolD_out=`otool -D #{dylibPath}`.strip.split(/\n/)
     has_install_name = otoolD_out.length() > 1
     otool_out=`otool -L #{dylibPath}`.strip.split(/\n/)

     # update install_name (-id) of current DyLib
     if has_install_name
       # strips first two lines
       id, otool_out = extractInstallName(otool_out)
       fixId(newDyLibPath, libname)
       debug "Fix install_name of DYLIB #{dylibPath} --> #{id}"
     else
       otool_out.delete_at(0).strip
     end

     # check the actual dependencies (lines 3++)
     handleDependencies(otool_out, targetPath, newDyLibPath, rpaths)

     # mark as processed
     $handledLibraries.add(libname)
   else
     debug "Already fixed #{libname}. Only change load path."
   end

   # readjust
   $currentIndent-=1

   return newDyLibPath,libname

  else

   debug "RPath but no lib found. Assuming fixed already: #{dylibPath}"
   libname=File.basename(dylibPath)
   $currentIndent-=1

   newDyLibPath="#{targetPath}/#{libname}"
   return newDyLibPath,libname

  end
end

###############################################################################
def handleBinary(binaryPath)
  $currentIndent+=1

  debug "Fixing binary #{binaryPath}"

  rpath = %x[otool -l #{binaryPath} | grep -A2 LC_RPATH | grep path | sed -n "s/^.*path\\s*\\(.*\\)(offset.*\$/\\1/p"]
  rpaths = rpath.split(/\n/)

  # no copy, no id change; juts run otool
  otool_out=`otool -L #{binaryPath}`.strip.split(/\n/)

  # remove 1st line containing the binary name
  otool_out.delete_at(0)

  # handle all referenced libraries and copy them if necessary to
  # the lib ref path
  # TODO parse RPATH?
  handleDependencies(otool_out, $lib_dir, binaryPath, rpaths)

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
  [ '--lib-path', '-l', GetoptLong::OPTIONAL_ARGUMENT ],
  [ '--install-name-tool', '-i', GetoptLong::OPTIONAL_ARGUMENT ],
  [ '--bin-path', '-b', GetoptLong::OPTIONAL_ARGUMENT ],
  [ '--plugin-path', '-p', GetoptLong::OPTIONAL_ARGUMENT ],
  [ '--path-prefix', '-e', GetoptLong::OPTIONAL_ARGUMENT ],
  [ '--extract-from-framework', '-f', GetoptLong::NO_ARGUMENT],
  [ '--no-auto-relative', '-n', GetoptLong::NO_ARGUMENT],
  [ '--no-copy', '-c', GetoptLong::NO_ARGUMENT]
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
  the path were optional plugin libraries like from Qt are located (experimental)
-e, --path-prefix:
  the prefix that is added to the new install_name (default: @executable_path/)
-f, --extract-from-framework:
  extract the linked libraries from their Framework folder. CAUTION: this does not copy Headers and Resources. Not tested with already present frameworks.
-n, --no-auto-relative:
  use path-prefix as is, without adding the relative path between bin path and lib path.
-c, --no-copy:
  do not copy dependencies, only fix up loading to relative paths
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
    when '--no-auto-relative'
      $auto_relative = false
    when '--no-copy'
      $nocopy = true
    when '--verbose'
      $DEBUG = true
  end
end

if $lib_dir.nil? and $bin_dir.nil?
  puts "Please provide at least a lib path or bin path with no-copy"
  puts usage.to_s
  exit 1
elsif $lib_dir.nil? and not $bin_dir.nil? and not $nocopy
  puts "If you only specify a bin_dir, no-copy must be active."
  puts usage.to_s
  exit 1
elsif !$lib_dir.nil? and !$bin_dir.nil?
  if $auto_relative
    $executableId = $executableId + $lib_dir.relative_path_from($bin_dir).to_s
  end
  if not $executableId.empty?
    $executableId += "/"
  end
  puts "Substituting prefix to find libs with:"
  puts $executableId
end


# fix libraries contained in lib-path
# recurse two-levels to capture the libraries inside the frameworks
if !$lib_dir.nil?
    debug "HANDLING LIB DIR"
    Dir.chdir($lib_dir.to_s) do
      lib_files = Dir.glob(["*","*/*"])
      for content in lib_files
        if ($plugin_dir.nil? or not ($lib_dir+content).to_s.start_with?($plugin_dir.to_s))
          if fixable(content, $lib_dir)
              if isFramework(content)
                handleFramework($lib_dir + content, $lib_dir, [])
              elsif (content.end_with?(".dylib") or content.end_with?(".so"))
                # TODO what to do with extracted Qt libs. They dont have an ending.
                handleDyLib($lib_dir + content, $lib_dir, [])
              else
                debug "Skipped #{$lib_dir + content} -- No lib or framework?"
              end
            else
              debug "Skipped #{$lib_dir + content} -- Otool not executable on it."
            end
        else
          debug "Skipped #{$lib_dir + content} -- Part of plugin_dir."
        end
      end
    end
end

# fix binary references
if !$bin_dir.nil?
  debug "HANDLING BIN DIR"
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

## TODO does not really work yet. Would copy all plugin dependencies into plugin_dir
##  but it SHOULD copy them into lib and relink relatively if necessary.
##  Currently we just copy the plugins and hope they have dependencies that are covered
##  by the usual libs. Should work for Qt (our only plugin).
if !$plugin_dir.nil?
  debug "HANDLING PLUGIN DIR"
  for content in Dir.glob("#{$plugin_dir}/**/*.dylib")
    debug "Handle dylib #{$plugin_dir + content}"
    handleDyLib($plugin_dir + content, File.dirname($plugin_dir + content), [])
  end
end


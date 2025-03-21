#!/usr/bin/ruby

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

require 'getoptlong'
require 'pathname'
require 'set'
require 'fileutils'

##### GLOBAL VARIABLES
$dir = nil
$codesign_tool=`which codesign`.strip
$signing_id = nil
$DEBUG = false

###############################################################################
def debug(message)
  if($DEBUG)
    puts "#{message}"
  end
end

###############################################################################
def deletePRL(path)
  if path.to_s.match(/.*\.app\/.*/)
    return
  end
  if path.to_s.match(/.*\.framework\/.*\.prl/)
    puts "Deleting #{path}"
    File.delete(path)
  end
end

def fixable(path)
  if path.to_s.match(/.*\.app\/.*/)
    return false
  elsif path.to_s.match(/.*\.app$/)
    return false
  elsif path.to_s.match(/.*\.framework$/)
    return true
  elsif path.to_s.match(/.*\.framework\/.*/)
    return false
  elsif path.to_s.match(/.*\.dylib$/)
    return true
  elsif path.to_s.match(/.*\.so$/)
    return true
  elsif File.directory?(path)
    return false
  elsif File.executable?(path)
    return true
  #elsif path.to_s.match(/.*share\/doc\/.*/)
  #  return false
  else
    return false
  end
end

def sign(path)
    `#{$codesign_tool} --deep --force --options runtime -s "#{$signing_id}" "#{path}"`
end

###############################################################################
###############################################################################
###############################################################################

##### INSTANTIATE CMD PARSER
opts = GetoptLong.new(
  [ '--help', '-h', GetoptLong::NO_ARGUMENT ],
  [ '--verbose', '-v', GetoptLong::NO_ARGUMENT ],
  [ '--dir-path', '-d', GetoptLong::REQUIRED_ARGUMENT ],
  [ '--signing-identity', '-s', GetoptLong::REQUIRED_ARGUMENT ],
  [ '--codesign-tool', '-c', GetoptLong::OPTIONAL_ARGUMENT ]
)

usage = "#{File.basename($0)} --bin-path PATH-TO-BINARIES --lib-path PATH-TO-STORE-LIBRARIES
-h, --help:
  show help
-c, --codesign-tool:
  specify the codesign-tool manually
-s, --signing-identity:
  specify the signing identity. Needs to be in default keychain and unlocked.
-d, --dir-path:
  the target path where the files should be signed
-v, --verbose:
  increase verbosity
"

opts.each do |opt, arg|
  case opt
    when '--help'
      puts usage.to_s
    when '--codesign-tool'
      $codesign_tool = arg
    when '--dir-path'
      $dir = Pathname.new(arg).realpath
    when '--signing-identity'
      $signing_id = arg
    when '--verbose'
      $DEBUG = true
  end
end

if $dir.nil?
  puts "Please provide a path"
  puts usage.to_s
  exit 1
end

if $signing_id.nil?
  puts "Please provide a signing identity"
  puts usage.to_s
  exit 1
end

debug "CLEANING frameworks"
for content in Dir.glob("#{$dir}/**/*")
  deletePRL($dir + content)
end

debug "HANDLING DIR"
for content in Dir.glob("#{$dir}/**/*")
  if fixable($dir + content)
    sign($dir + content)
  end
end

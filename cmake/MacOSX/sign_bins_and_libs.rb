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

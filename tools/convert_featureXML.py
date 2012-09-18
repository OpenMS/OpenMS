#!/usr/bin/env python
# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry               
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.

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
# $Maintainer: Hendrik Weisser $
# $Authors: Hendrik Weisser $
# --------------------------------------------------------------------------

"""Convert featureXML files between formats used before and after OpenMS 1.8"""

import sys
import re
import os.path

if len(sys.argv) < 3:
    sys.exit("Error: two parameters (input featureXML, output featureXML) "
             "expected!")

in_path = sys.argv[1]
out_path = sys.argv[2]

if os.path.realpath(in_path) == os.path.realpath(out_path):
    sys.exit("Error: output path must be different from input path!")

old_pattern = re.compile('<hullpoint>\n<hposition dim="0">(.*)</hposition>\n'
                         '<hposition dim="1">(.*)</hposition>\n</hullpoint>')
new_pattern = re.compile('<pt x="(.*)" y="(.*)" */>')

with open(in_path) as source:
    with open(out_path, "w") as sink:
        line_count = 0
        for line in source:
            line_count += 1
            stripped = line.strip()
            if stripped == "<hullpoint>": # convert old to new:
                temp_count = line_count
                indent = line[:line.index("<")]
                content = stripped
                while stripped != "</hullpoint>":
                    stripped = source.next().strip()
                    line_count += 1
                    content += "\n" + stripped
                match = old_pattern.match(content)
                if not match:
                    print "skipping unexpected content in lines %d-%d" % \
                        (temp_count, line_count)
                else:
                    sink.write(indent + '<pt x="%s" y="%s" />\n' % 
                               (match.group(1), match.group(2)))
                continue
            # else:
            match = new_pattern.match(stripped)
            if match: # convert new to old:
                indent = line[:line.index("<")]
                sink.write(indent + "<hullpoint>\n")
                sink.write(indent + '\t<hposition dim="0">%s</hposition>\n' %
                           match.group(1))
                sink.write(indent + '\t<hposition dim="1">%s</hposition>\n' %
                           match.group(2))
                sink.write(indent + "</hullpoint>\n")
            else:
                sink.write(line)

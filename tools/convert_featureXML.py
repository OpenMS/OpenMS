#!/usr/bin/env python

"""Convert featureXML files between formats used before and after OpenMS 1.8"""

import sys
import re

if len(sys.argv) < 3:
    sys.exit("two parameters (input featureXML, output featureXML) expected!")

in_path = sys.argv[1]
out_path = sys.argv[2]

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

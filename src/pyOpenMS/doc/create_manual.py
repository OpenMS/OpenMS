import os
import glob
import re
pxd_files = glob.glob("../pxds/*.pxd")
classdocu_base = "http://ftp.mi.fu-berlin.de/OpenMS/release-documentation/html"
classdocu_base = "http://www.openms.de/current_doxygen/html/"

# Variables:
# after how many characters should a line be wrapped
Linebreak = 85
# Whether to remove comments and "nogil" statements
remove_nogil = True
remove_comments = True
ignores = ["ctime", "smart_ptr", "Time"]

# Regexes to determine if a line is comment only or empty
comment_line = re.compile("^\s*#")
empty_line = re.compile("^\s*$")

# Some state variables to keep track of 
#  - whether the previous line was empty or not
#  - whether the current code block contains any pyOpenMS specific "wrap-XXX" statements
# 
# We want to remove all comments except those pertaining to autowrap (Python
# wrapping) directives as they are important to understand the code
prevline_empty = False
wrap_found = False

# Apply line break
def line_break(line):
    # print "trying to break line ", line
    if len(line) < Linebreak:return line

    for i in reversed(range(Linebreak)):
        if line[i] == " ":
            break
    if len(line) - i < 9: return line
    if i == 0:
        print "Could not break line", line
        # return line
        return line[:Linebreak] + "\\\n" + " "*6 + line_break(line[Linebreak:])
    # print "call again with i", i, "and len ", len(line)
    return line[:i] + "\\\n" + " "*6 + line_break(line[i:])

# Prepare a line for printing
def prepare_line(line):
    global wrap_found, prevline_empty
    # replace the indent from 4/8 to 2/4 -> get more space
    line = line.replace("        ", "xxxx")
    line = line.replace("    ", "  ")
    line = line.replace("xxxx", "    ")
    if remove_nogil:
        line = line.replace("nogil except +", "")

    # reset
    if prevline_empty:
        wrap_found = False

    if remove_comments:
        if line.find("wrap-") != -1:
            wrap_found = True
        if re.search(comment_line, line) is not None and not wrap_found:
            if line.find("COMMENT:") != -1:
                line = line.replace("COMMENT:", "")
            elif line.find(" COMMENT:") != -1:
                line = line.replace(" COMMENT:", "")
            else:
                print "Remove line", line.strip()
                return ""
    return line

def get_namespace(pxd):
    filehandle = open(pxd)
    fulltext = filehandle.read()
    filehandle.close()
    import re
    match = re.search("cdef extern.*?namespace\s*\"([^\"]*)\"", fulltext)
    if not match:
        return "OpenMS"
    else:
        return match.group(1)


def get_header(title, ns):

    return r"""
$\rightarrow$ \textit{\href{%s/class%s_1_1%s.html}{Link
to OpenMS documentation}}

Wrapped functions in Python:
""" % (classdocu_base, ns, title)

write_handle = open("appendix.tex", "w")
for pxd in sorted(pxd_files):

    # get title for latex section
    wrap_found = False
    title = pxd[8:-4]
    if title in ignores:
        continue

    title = title.replace("_", "\\_")
    write_handle.write("\subsection{%s}\n\n" % title)
    write_handle.write("\label{%s}\n\n" % title)
    # write_handle.write("{\\tiny\n    \\begin{verbatim}")

    # Write the header (link to OpenMS documentation)
    write_handle.write( get_header(title, get_namespace(pxd)) )
    write_handle.write("{    \\begin{verbatim}")

    # Start writing the PDX file content
    filehandle = open(pxd)
    found_start = False
    for line in filehandle:
        if line.find("import") == -1 and len(line.strip() ) > 0:
            found_start = True
        if line.find("import") == -1 and found_start:
            # print line.strip()
            line = prepare_line(line)
            line = line_break(line)

            # do not print two empty lines in a row
            if re.search(empty_line, line) is not None:
                if prevline_empty:
                    continue
                prevline_empty = True
            else:
                prevline_empty = False

            write_handle.write(line)

    write_handle.write("\end{verbatim}\n}\n\n")

write_handle.close()




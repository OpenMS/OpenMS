import os
import glob
pxd_files = glob.glob("../pxds/*.pxd")

# after how many characters should a line be wrapped
Linebreak = 85
remove_nogil = True
ignores = ["ctime", "smart_ptr", "Time"]

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
    # replace the indent from 4/8 to 2/4 -> get more space
    line = line.replace("        ", "xxxx")
    line = line.replace("    ", "  ")
    line = line.replace("xxxx", "    ")
    return line

write_handle = open("appendix.tex", "w")
print pxd_files
for pxd in sorted(pxd_files):
    # get title for latex section
    title = pxd[8:-4]
    if title in ignores:
        continue
    title = title.replace("_", "\\_")
    write_handle.write("\subsection{%s}\n\n" % title)
    write_handle.write("\label{%s}\n\n" % title)
    # write_handle.write("{\\tiny\n    \\begin{verbatim}")
    write_handle.write("{    \\begin{verbatim}")
    filehandle = open(pxd)
    found_start = False
    for line in filehandle:
        if line.find("import") == -1 and len(line.strip() ) > 0:
            found_start = True
        if line.find("import") == -1 and found_start:
            # print line.strip()
            line = prepare_line(line)
            line = line_break(line)
            write_handle.write(line)
    write_handle.write("\end{verbatim}\n}\n\n")

write_handle.close()




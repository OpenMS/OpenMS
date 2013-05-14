import os
import glob
pxd_files = glob.glob("../pxds/*.pxd")

# after how many characters should a line be wrapped
Linebreak = 85
remove_nogil = True

# Apply line break
def line_break(line):
    if len(line) < Linebreak:return line

    for i in reversed(range(Linebreak)):
        if line[i] == " ":
            break
    if len(line) - i < 9: return line
    return line[:i] + "\\\n" + " "*6 + line_break(line[i:])

# Prepare a line for printing
def prepare_line(line):
    # replace the indent from 4/8 to  2/4 -> get more space
    line = line.replace("        ", "xxxx")
    line = line.replace("    ", "  ")
    line = line.replace("xxxx", "    ")
    return line

write_handle = open("appendix.tex", "w")
for pxd in pxd_files:
    pxd
    title = pxd[8:-4]
    title = title.replace("_", "\\_")
    write_handle.write("\subsection{%s}\n\n" % title)
    # write_handle.write("{\\tiny\n    \\begin{verbatim}")
    write_handle.write("{    \\begin{verbatim}")
    filehandle = open(pxd)
    for line in filehandle:
        if line.find("import") == -1:
            print line.strip()
            line = prepare_line(line)
            line = line_break(line)
            write_handle.write(line)
    write_handle.write("\end{verbatim}\n}\n\n")

write_handle.close()




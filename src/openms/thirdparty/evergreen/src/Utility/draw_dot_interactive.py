import pygraphviz as pgv
import sys
import pylab as P
#P.ion()

def main(argv):
    if len(argv) != 2:
        print 'Usage: draw_dot.py <graph .dot> <graph .png output>'
    else:
        G = pgv.AGraph(argv[0])

        G.layout(prog='circo')
        G.layout(prog='neato', args='-Goverlap=scale -Gsplines=true')
        G.draw(argv[1])

if __name__ == "__main__":
    main(sys.argv[1:])

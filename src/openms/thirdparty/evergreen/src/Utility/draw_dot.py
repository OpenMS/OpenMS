import pygraphviz as pgv
import sys

def main(argv):
    if len(argv) not in (2,3,4,5):
        print 'Usage: draw_dot.py <graph .dot> <output image file> [layout, default is neato; optionally multiple layouts separated by ","] [overlap mode {scale, false}, default is false] [label edges {0,1}]'
    else:
        layouts = ['neato']
        if len(argv) >= 3:
            layouts = argv[2].split(',')

        G = pgv.AGraph(argv[0])

        overlap='false'
        if len(argv) >= 4:
            overlap=argv[3]

        label_edges=True
        if len(argv) >= 5:
            label_edges = bool(int(argv[4]))
        if not label_edges:
            for e in G.edges():
                e.attr['label'] = ' '
        
        for lo in layouts[:-1]:
            G.layout(prog=lo)
        G.layout(prog=layouts[-1], args='-Goverlap=' + overlap + ' -Gsplines=true')

        G.draw(argv[1])

if __name__ == "__main__":
    main(sys.argv[1:])

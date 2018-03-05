from __future__ import print_function

def parse_pxd_file(path):
    import pdb

    import os

    from Cython.Compiler.CmdLine import parse_command_line
    from Cython.Compiler.Main import create_default_resultobj, CompilationSource
    from Cython.Compiler import Pipeline
    from Cython.Compiler.Scanning import FileSourceDescriptor
    from Cython.Compiler.Nodes import CEnumDefNode, CppClassNode, CTypeDefNode, CVarDefNode, CImportStatNode, CDefExternNode
    # from Cython.Compiler.ExprNodes import *


    options, sources = parse_command_line(["--cplus", path])

    path = os.path.abspath(path)
    basename = os.path.basename(path)
    name, ext = os.path.splitext(basename)

    source_desc = FileSourceDescriptor(path, basename)
    source = CompilationSource(source_desc, name, os.getcwd())
    result = create_default_resultobj(source, options)

    context = options.create_context()
    pipeline = Pipeline.create_pyx_pipeline(context, options, result)
    context.setup_errors(options, result)
    root = pipeline[0](source)  # only parser

    def iter_bodies(tree):
        try:
            for n in tree.body.stats[0].stats:
                # cimports at head of file
                yield n
        except:
            pass
        if hasattr(tree.body, "stats"):
            for s in tree.body.stats:
                if isinstance(s, CDefExternNode):
                    body = s.body
                    if hasattr(body, "stats"):
                        for node in body.stats:
                            yield node
                    else:
                        yield body
        elif hasattr(tree.body, "body"):
            body = tree.body.body
            yield body
        else:
            raise Exception("parse_pxd_file failed: no valied .pxd file !")

    lines = open(path).readlines()

    def cimport(b, _, __):
        print ("cimport", b.module_name, "as", b.as_name)

    handlers = { CEnumDefNode : "",
                 CppClassNode : "",
                 CTypeDefNode : "",
                 CVarDefNode  : "",
                 CImportStatNode  : "",
                 }

    # from autowrap.PXDParser import EnumDecl.parseTree, CppClassDecl.parseTree, CTypeDefDecl.parseTree, MethodOrAttributeDecl.parseTree
    # handlers = { CEnumDefNode : EnumDecl.parseTree,
    #              CppClassNode : CppClassDecl.parseTree,
    #              CTypeDefNode : CTypeDefDecl.parseTree,
    #              CVarDefNode  : MethodOrAttributeDecl.parseTree,
    #              CImportStatNode  : cimport,
    #              }

    result = []
    for body in iter_bodies(root):
        handler = handlers.get(type(body))
        if handler is not None:
            # L.info("parsed %s, handler=%s" % (body.__class__, handler.im_self))
            result.append([ body, lines, path ])
        else:
            for node in getattr(body, "stats", []):
                handler = handlers.get(type(node))
                if handler is not None:
                    result.append([ node, lines, path ])
    return result


def create_pxd_file_map(src_path):
    # Searches through all pxd files and matches them to the corresponding
    # OpenMS headers. Note that each header can be mapped to multiple pxd
    # files:
    #
    # pxd_file_matching[ HeaderName ] = set( [PXDFile, PXDFile, ...] )
    #
    import glob, os, re
    pxd_path = os.path.join(src_path, "src/pyOpenMS/pxds/")
    pxds = glob.glob(pxd_path + "/*.pxd")
    pxd_file_matching = {}
    for pfile in pxds:
        filename_rgx = re.compile("cdef extern from \"<([^\"]*)>\"", re.IGNORECASE)
        filematch = [o.group(1) for o in filename_rgx.finditer(open(pfile).read()) ]
        for fm in filematch:
            res = pxd_file_matching.get(fm, set([]))
            res.update([pfile])
            pxd_file_matching[fm] = res

    return pxd_file_matching


from docutils import nodes
from docutils.parsers.rst.roles import set_classes

def path_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    """Create path

    Returns 2 part tuple containing list of nodes to insert into the
    document and a list of system messages.  Both are allowed to be
    empty.

    :param name: The role name used in the document.
    :param rawtext: The entire markup snippet, with role.
    :param text: The text marked with the role.
    :param lineno: The line number where rawtext appears in the input.
    :param inliner: The inliner instance that called us.
    :param options: Directive options for customization.
    :param content: The directive content for customization.
    """
    try:
        pathparts = text.split(",")
        if len(pathparts) == 0:
            raise ValueError
    except ValueError:
        msg = inliner.reporter.error(
            'Path length must be greater or equal to 1; '
            '"%s" is invalid.' % text, line=lineno)
        prb = inliner.problematic(rawtext, rawtext, msg)
        return [prb], [msg]
    app = inliner.document.settings.env.app
    nodes = make_path_nodes(rawtext, pathparts, app, options)
    return nodes, []

def make_path_nodes(rawtext, pathparts, app, options):
    """Create a path with folder image and separators

    :param rawtext: Text being replaced with link node.
    :param app: Sphinx application context
    :param pathparts: The parts of the path
    :param options: Options dictionary passed to role func.
    """
    #
    try:
        base = app.config.pathicon
        if not base:
            raise AttributeError
    except AttributeError:
        raise ValueError('pathicon configuration value is not set')
    #
    triangle = nodes.inline(text="▸")
    #triangle = nodes.Text("▹")
    triangle["classes"].append("pathsep_triangle")

    #retnodes = [nodes.image(app.config.pathicon, **options), triangle]
    nobrhtml = "<nobr class='path'>"
    nobrelem = nodes.raw('', nobrhtml, format="html")
    brhtml = "<br>"
    brelem = nodes.raw('', brhtml, format="html")
    nobrhtmlclose = "</nobr>"
    nobrelemclose = nodes.raw('', nobrhtmlclose, format="html")

    html = f"<i class=\"{app.config.pathicon} path\"></i>"
    retnodes = [nobrelem, nodes.raw('', html, format="html")]
    for part in pathparts[:-1]:
        if part == "/break":
            retnodes.append(nobrelemclose)
            retnodes.append(brelem)
            retnodes.append(nobrelem)
            retnodes.append(triangle.deepcopy())
        else:
            retnodes.append(nodes.Text(part))
            retnodes.append(triangle.deepcopy())
    retnodes.append(nodes.Text(pathparts[-1]))
    retnodes.append(nobrelemclose)
    #slash = '/' if base[-1] != '/' else ''
    #ref = base + slash + type + '/' + slug + '/'
    set_classes(options)
    #node = nodes.reference(rawtext, type + ' ' + utils.unescape(slug), refuri=ref,
    #                       **options)
    return retnodes

def setup(app):
    """Install the plugin.

    :param app: Sphinx application context.
    """
    app.add_role('path', path_role)
    app.add_config_value('pathicon', None, 'env')
    return

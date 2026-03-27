from docutils.nodes import math


def chem_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    latex = rf'\ce{{{text}}}'
    node = math(rawtext, latex, **options)
    return [node], []


def setup(app):
    """Install the plugin.
    :param app: Sphinx application context.
    """
    app.add_role('chem', chem_role)
    return {'parallel_read_safe': True}

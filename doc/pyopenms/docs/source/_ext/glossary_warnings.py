from docutils import nodes
from typing import cast
from sphinx import addnodes
from sphinx.util import logging

logger = logging.getLogger(__name__)


class FindTextNodesVisitor(nodes.SparseNodeVisitor):

    def __init__(self, document, words):
        super().__init__(document)
        self.words = words

    def visit_Text(self, node):
        for substring in self.words:
            if substring in node.astext():
                logger.warn(
                    logging.get_node_location(node) +
                    ' Potential glossary terms found in the text node. First match:' +
                    substring)
                break

    def visit_literal_block(self, node):
        self.in_literal_block = True
        raise nodes.SkipChildren

    def visit_term(self, node):
        self.in_term = True
        raise nodes.SkipChildren


def collect_glossary_terms(app, doctree):
    env = app.builder.env

    if not hasattr(env, 'glossary_all_terms'):
        env.glossary_all_terms = set()
    for glossary in doctree.findall(addnodes.glossary):
        definition_list = cast(nodes.definition_list, glossary[0])
        for t in definition_list:
            env.glossary_all_terms.update(
                [" " + w + " "
                    for w in cast(nodes.term, t).astext().split("\n\n")[0].split("\n")])


def check_forbidden_words(app, doctree, docname):
    if ("_autosummary" not in docname):
        env = app.builder.env
        if not hasattr(env, 'glossary_all_terms'):
            env.glossary_all_terms = set()
        visitor = FindTextNodesVisitor(doctree, env.glossary_all_terms)
        doctree.walk(visitor)


def purge_glossary_terms(app, env, docname):
    if not hasattr(env, 'glossary_all_terms'):
        return

    env.glossary_all_terms = [term for term in env.glossary_all_terms
                              if term['docname'] != docname]


def merge_glossary_terms(app, env, docnames, other):
    if not hasattr(env, 'glossary_all_terms'):
        env.glossary_all_terms = set()
    if hasattr(other, 'glossary_all_terms'):
        env.glossary_all_terms.update(other.glossary_all_terms)


def setup(app):
    app.connect('doctree-resolved', check_forbidden_words)
    app.connect('doctree-read', collect_glossary_terms)
    app.connect('env-merge-info', merge_glossary_terms)
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------

project = 'OpenMS'
copyright = '2002-present, OpenMS Inc'
author = 'OpenMS Team'


# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short major.minor.patch version.
version = '3.5.0'
# Short version for the latest supported KNIME
knime_version = '5.5.0'

# The full version, including alpha/beta/rc tags.
release = '3.5.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
import sys, os

sys.path.append(os.path.abspath('../_ext'))

extensions = [
  'sphinx_copybutton',
  'sphinx.ext.autodoc',
  'sphinx.ext.autosummary',
  'myst_parser',
  'notfound.extension',
  'sphinxcontrib.images',
  'sphinx_inline_tabs',
  'pathrole',
  'hoverxref.extension',
  'sphinx_search.extension',
  'sphinx_design',
]

numfig = True

myst_enable_extensions = [
  "tasklist",
  "dollarmath",
  "amsmath",
  "colon_fence",
  "linkify",
  "replacements",
  "html_admonition",
  "substitution",
]

root_doc = 'index'
myst_heading_anchors = 3

# do not start footnotes with transition, gives warnings/errors when
# foot notes are directly used after a heading.
myst_footnote_transition = False

autosummary_generate = True
autosummary_imported_members = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = ['.rst', '.md']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store',
    'openms-applications-and-tools/installation/installation-with-conda.md', # just a snippet to be included
    'openms-applications-and-tools/installation/run-in-container.md', # just a snippet to be included 
]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'furo'

html_static_path = ['../_static']
html_favicon = '../assets/logo/favicon.png'
html_theme_options = {
    "navigation_with_keys": True,
    "light_css_variables": {
        "font-stack--monospace": "Consolas, monospace",
        "font-size--small": "90%",
        "toc-font-size": "87.5%"
    },
    "light_logo": "FinalLogo_Nov2024_Versions-01.svg",
    "dark_logo": "FinalLogo_Nov2024_Versions-02.svg",
}
pygments_style = 'sas'

pygments_dark_style = 'rrt'

# Configure tooltips
hoverxref_roles = ['term']

hoverxref_role_types = {'term':'tooltip'}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

html_css_files = [
    'css/custom.css',
    'css/all.min.css',
]

html_js_files = [
    'js/matomo.js',  # Matomo tracking code for openms.matomo.cloud
    'js/piwik.js',  # For tracking user statistics on openms-web.piwik.pro
    'js/gmc.min.js',  # GitHub Manyfaced Cards for showcasing GitHub repos (e.g. nf-core workflows)
]

#pathicon = 'data:image/png;base64,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'

pathicon = 'fa fa-folder-open'

variables_to_export = [
    "project",
    "version",
    "knime_version"
]
myst_substitutions = {}
for v in variables_to_export:
    myst_substitutions[v] = globals()[v]

frozen_locals = dict(locals())
# Makes variables_to_export available in the epilog
rst_epilog = '\n'.join(map(lambda x: f".. |{x}| replace:: {frozen_locals[x]}", variables_to_export))
del frozen_locals

"""Make Sphinx autodoc/autosummary recognise nanobind functions.

pyOpenMS is built with nanobind, which binds methods as ``nanobind.nb_method``
and static methods and free functions as ``nanobind.nb_func``. Neither type is
a ``builtin_function_or_method``, so autodoc misclassifies them:

* ``nb_method`` is a non-data descriptor, which autodoc's AttributeDocumenter
  claims with a higher priority than MethodDocumenter -- every method ends up
  in the "Attributes" table of a class page.
* ``nb_func`` is not even a routine, so it is documented as an attribute whose
  value is rendered as ``<nanobind.nb_func object>``.

This extension routes both types to the method/function documenters. The
signatures come from the first lines of the nanobind docstring (enabled via
``autodoc_docstring_signature``).
"""

import re

from sphinx.ext.autodoc import (
    AttributeDocumenter,
    ClassDocumenter,
    DataDocumenter,
    FunctionDocumenter,
    MethodDocumenter,
    ModuleDocumenter,
)

_NB_TYPES = ("nb_func", "nb_method")

# nanobind names bound types after the extension module that defines them
# (pyopenms._pyopenms_kernel.MSSpectrum); they are documented as pyopenms.X
_PRIVATE_MODULE = re.compile(r"\bpyopenms\._pyopenms_\w+\.")
# leading ``self`` of an instance method signature
_SELF_ARG = re.compile(r"^\(self(?:\s*,\s*|(?=\)))", re.MULTILINE)


def is_nanobind_function(obj):
    tp = type(obj)
    return tp.__module__ == "nanobind" and tp.__name__ in _NB_TYPES


def is_nanobind_static(obj):
    # nb_method is reserved for functions bound with an implicit ``self``
    return is_nanobind_function(obj) and type(obj).__name__ == "nb_func"


def _bound_function(obj):
    """The nanobind function behind ``obj``, also through Python wrappers.

    pyopenms.addons wraps some bindings with ``functools.update_wrapper``;
    those wrappers carry the nanobind docstring (and its signatures).
    """
    for _ in range(5):
        if is_nanobind_function(obj):
            return obj
        obj = getattr(obj, "__wrapped__", None)
        if obj is None:
            break
    return None


def _clean_signature(sig, strip_self, summary):
    sig = _PRIVATE_MODULE.sub("pyopenms.", sig)
    if strip_self:
        sig = _SELF_ARG.sub("(", sig)
    if summary:
        # overloads yield one signature per line, which autosummary's
        # signature mangling cannot handle; the table shows the first one
        sig = sig.split("\n", 1)[0]
    return sig


def _is_summary(kwargs):
    # autosummary asks for signatures without annotations, autodoc does not
    return kwargs.get("show_annotation") is False


class NanobindMethodDocumenter(MethodDocumenter):
    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        if is_nanobind_function(member):
            return isinstance(parent, ClassDocumenter)
        return super().can_document_member(member, membername, isattr, parent)

    def import_object(self, raiseerror=False):
        ret = super().import_object(raiseerror)
        if ret and is_nanobind_static(self.object):
            # document static members before regular methods, like autodoc
            # does for Python staticmethods
            self.member_order = MethodDocumenter.member_order - 1
        return ret

    def add_directive_header(self, sig):
        super().add_directive_header(sig)
        if is_nanobind_static(self.object):
            self.add_line("   :staticmethod:", self.get_sourcename())

    def format_signature(self, **kwargs):
        sig = super().format_signature(**kwargs)
        target = _bound_function(self.object)
        if target is None:
            return sig
        return _clean_signature(
            sig,
            strip_self=not is_nanobind_static(target),
            summary=_is_summary(kwargs),
        )


class NanobindFunctionDocumenter(FunctionDocumenter):
    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        if is_nanobind_function(member):
            return isinstance(parent, ModuleDocumenter)
        return super().can_document_member(member, membername, isattr, parent)

    def format_signature(self, **kwargs):
        sig = super().format_signature(**kwargs)
        if _bound_function(self.object) is None:
            return sig
        return _clean_signature(sig, strip_self=False, summary=_is_summary(kwargs))


class NanobindAttributeDocumenter(AttributeDocumenter):
    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # also covers Python wrappers: autodoc unwraps them to the nanobind
        # function and would then take them for an attribute descriptor
        if _bound_function(member) is not None:
            return False
        return super().can_document_member(member, membername, isattr, parent)


class NanobindDataDocumenter(DataDocumenter):
    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        if _bound_function(member) is not None:
            return False
        return super().can_document_member(member, membername, isattr, parent)


def setup(app):
    app.setup_extension("sphinx.ext.autodoc")
    # same objtype as the stock documenters, so autosummary's "methods",
    # "functions" and "attributes" tables pick up the corrected classification
    for documenter in (
        NanobindMethodDocumenter,
        NanobindFunctionDocumenter,
        NanobindAttributeDocumenter,
        NanobindDataDocumenter,
    ):
        app.add_autodocumenter(documenter, override=True)
    return {"parallel_read_safe": True, "parallel_write_safe": True}

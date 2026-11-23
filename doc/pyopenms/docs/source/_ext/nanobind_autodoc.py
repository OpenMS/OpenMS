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

from docutils import nodes
from sphinx.ext.autodoc import (
    AttributeDocumenter,
    ClassDocumenter,
    DataDocumenter,
    FunctionDocumenter,
    MethodDocumenter,
    ModuleDocumenter,
    PropertyDocumenter,
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


def _getter_signature(prop):
    """``(self) -> T`` of a property bound by nanobind (def_rw, def_prop_rw)."""
    getter = getattr(prop, "fget", None)
    if not is_nanobind_function(getter):
        return None
    sigs = getter.__nb_signature__
    if len(sigs) != 1 or not sigs[0][0].startswith("def (self) -> "):
        return None
    return sigs[0][0][len("def ") :]


class NanobindPropertyDocumenter(PropertyDocumenter):
    def add_directive_header(self, sig):
        start = len(self.directive.result)
        super().add_directive_header(sig)
        getter_sig = _getter_signature(self.object)
        if getter_sig is None or self.config.autodoc_typehints == "none":
            return
        # autodoc cannot introspect the nanobind getter, so it adds no type
        if not any(
            line.lstrip().startswith(":type:") for line in self.directive.result[start:]
        ):
            objtype = getter_sig[len("(self) -> ") :]
            self.add_line("   :type: " + _public_name(objtype), self.get_sourcename())


class NanobindClassDocumenter(ClassDocumenter):
    def get_canonical_fullname(self):
        # the canonical name of a nanobind class is in a private extension
        # module (pyopenms._pyopenms_kernel.X); emitting it registered every
        # enum that is both nested and exported at module level twice
        canonical = super().get_canonical_fullname()
        if canonical and _PRIVATE_MODULE.match(canonical):
            return None
        return canonical

    def filter_members(self, members, want_all):
        # A nested class (mostly an enum such as ProgressLogger.LogType) is listed
        # by every subclass (:inherited-members:), so each copy documented the
        # same canonical name again. Keep it only on its defining class, which
        # is the name the nanobind signatures refer to.
        owner = getattr(self.object, "__qualname__", "")
        kept = []
        for name, member, isattr in super().filter_members(members, want_all):
            if isinstance(member, type):
                defined_in = getattr(member, "__qualname__", "").rpartition(".")[0]
                if defined_in and defined_in != owner:
                    continue
            kept.append((name, member, isattr))
        return kept


def _unwrap_singletons():
    """Document the DB singletons (ModificationsDB, ElementDB, ...) as classes.

    pyopenms replaces them by callable wrappers (``ModificationsDB()`` returns
    the instance), which autosummary does not recognise as classes, so they had
    no API page. The wrappers carry the class as ``_singleton_class``.
    """
    try:
        import pyopenms
    except ImportError:
        return
    for name in dir(pyopenms):
        wrapped = getattr(getattr(pyopenms, name), "_singleton_class", None)
        if isinstance(wrapped, type) and wrapped.__name__ == name:
            setattr(pyopenms, name, wrapped)


def _drop_placeholder_class_signature(
    app, what, name, obj, options, signature, return_annotation
):
    # autodoc cannot introspect an overloaded nanobind constructor and falls
    # back to the meaningless "(*args, **kwargs)"; the overloads are listed
    # under __init__
    if (
        what == "class"
        and signature == "(*args, **kwargs)"
        and type(obj).__module__ == "nanobind"
    ):
        return "", return_annotation
    return None


# "1. ``name(self, a: pyopenms._pyopenms_kernel.X) -> None``" in the body
# nanobind writes for overloaded functions
_OVERLOAD_LINE = re.compile(r"^(\d+\. ``\w+\()self(?:, |(?=\)))")


def _clean_docstring(app, what, name, obj, options, lines):
    if isinstance(obj, property):
        # nanobind gives an undocumented property the getter signature as its
        # docstring; the type is shown as the property's :type: instead
        getter_sig = _getter_signature(obj)
        if getter_sig and "\n".join(lines).strip() == getter_sig:
            lines.clear()
        return
    target = _bound_function(obj)
    if target is None:
        return
    strip_self = not is_nanobind_static(target)
    for i, line in enumerate(lines):
        line = _PRIVATE_MODULE.sub("pyopenms.", line)
        if strip_self:
            line = _OVERLOAD_LINE.sub(r"\1", line)
        lines[i] = line


def _public_name(text):
    return _PRIVATE_MODULE.sub("pyopenms.", text)


def _resolve_private_reference(app, env, node, contnode):
    """Resolve pyopenms._pyopenms_kernel.X as its public name pyopenms.X.

    autodoc names a class by the module that defines it in "Bases:" and in
    the "alias of" line of an alias such as PeakMap; for nanobind that is a
    private extension module, which is not documented.
    """
    target = node.get("reftarget", "")
    if node.get("refdomain") != "py" or not _PRIVATE_MODULE.match(target):
        return None
    if isinstance(contnode, nodes.Text):
        contnode = nodes.Text(_public_name(contnode))
    else:
        for text in list(contnode.findall(nodes.Text)):
            text.parent.replace(text, nodes.Text(_public_name(text)))
    return env.get_domain("py").resolve_xref(
        env,
        node["refdoc"],
        app.builder,
        node["reftype"],
        _public_name(target),
        node,
        contnode,
    )


def setup(app):
    app.setup_extension("sphinx.ext.autodoc")
    _unwrap_singletons()
    app.connect("autodoc-process-signature", _drop_placeholder_class_signature)
    app.connect("autodoc-process-docstring", _clean_docstring)
    app.connect("missing-reference", _resolve_private_reference)
    # same objtype as the stock documenters, so autosummary's "methods",
    # "functions" and "attributes" tables pick up the corrected classification
    for documenter in (
        NanobindClassDocumenter,
        NanobindMethodDocumenter,
        NanobindFunctionDocumenter,
        NanobindAttributeDocumenter,
        NanobindDataDocumenter,
        NanobindPropertyDocumenter,
    ):
        app.add_autodocumenter(documenter, override=True)
    return {"parallel_read_safe": True, "parallel_write_safe": True}

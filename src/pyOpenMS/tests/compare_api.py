#!/usr/bin/env python3
"""
Compare two pyopenms API JSON files (from extract_api.py) and produce a report.

Usage:
    python compare_api.py reference.json current.json [--output report.md]

The reference is typically the Cython 3.5 API, and current is the nanobind 3.6 API.
"""

import json
import sys
import argparse


def load_api(path):
    with open(path) as f:
        return json.load(f)


def compare_methods(ref_methods, cur_methods, class_name):
    """Compare methods between reference and current."""
    issues = []
    ref_names = set(ref_methods.keys())
    cur_names = set(cur_methods.keys())

    missing = sorted(ref_names - cur_names)
    added = sorted(cur_names - ref_names)

    if missing:
        issues.append(('missing_methods', missing))
    if added:
        issues.append(('added_methods', added))

    return issues


def compare_classes(ref_cls, cur_cls, class_name):
    """Compare two class definitions."""
    issues = []

    # Compare methods
    ref_methods = ref_cls.get('methods', {})
    cur_methods = cur_cls.get('methods', {})
    method_issues = compare_methods(ref_methods, cur_methods, class_name)
    issues.extend(method_issues)

    # Compare properties
    ref_props = set(ref_cls.get('properties', []))
    cur_props = set(cur_cls.get('properties', []))
    missing_props = sorted(ref_props - cur_props)
    added_props = sorted(cur_props - ref_props)
    if missing_props:
        issues.append(('missing_properties', missing_props))
    if added_props:
        issues.append(('added_properties', added_props))

    # Compare nested enums
    ref_enums = set(ref_cls.get('nested_enums', {}).keys())
    cur_enums = set(cur_cls.get('nested_enums', {}).keys())
    missing_enums = sorted(ref_enums - cur_enums)
    if missing_enums:
        issues.append(('missing_nested_enums', missing_enums))

    # Compare nested enum values
    for enum_name in sorted(ref_enums & cur_enums):
        ref_values = set(ref_cls['nested_enums'][enum_name].get('values', {}).keys())
        cur_values = set(cur_cls.get('nested_enums', {}).get(enum_name, {}).get('values', {}).keys())
        missing_vals = sorted(ref_values - cur_values)
        if missing_vals:
            issues.append(('missing_enum_values', {enum_name: missing_vals}))

    # Compare static methods
    ref_static = set(ref_cls.get('static_methods', []))
    cur_static = set(cur_cls.get('static_methods', []))
    # Methods that were static in ref but not in current (or vice versa)
    ref_all_methods = set(ref_cls.get('methods', {}).keys())
    cur_all_methods = set(cur_cls.get('methods', {}).keys())
    common_methods = ref_all_methods & cur_all_methods
    static_changes = []
    for m in sorted(common_methods):
        was_static = m in ref_static
        is_static = m in cur_static
        if was_static != is_static:
            static_changes.append(f"{m}: {'static' if was_static else 'instance'} -> {'static' if is_static else 'instance'}")
    if static_changes:
        issues.append(('static_changes', static_changes))

    return issues


def compare_enums(ref_enum, cur_enum, enum_name):
    """Compare two enum definitions."""
    issues = []
    ref_values = set(ref_enum.get('values', {}).keys())
    cur_values = set(cur_enum.get('values', {}).keys())

    missing = sorted(ref_values - cur_values)
    added = sorted(cur_values - ref_values)

    if missing:
        issues.append(('missing_values', missing))
    if added:
        issues.append(('added_values', added))

    # Check for value changes
    for name in sorted(ref_values & cur_values):
        ref_val = ref_enum['values'][name]
        cur_val = cur_enum['values'][name]
        if str(ref_val) != str(cur_val):
            issues.append(('value_changed', f"{name}: {ref_val} -> {cur_val}"))

    return issues


def generate_report(ref, cur):
    """Generate a comparison report."""
    lines = []

    ref_version = ref['_meta']['version']
    cur_version = cur['_meta']['version']

    lines.append(f"# pyOpenMS API Comparison Report")
    lines.append(f"")
    lines.append(f"**Reference (Cython):** {ref_version}")
    lines.append(f"**Current (Nanobind):** {cur_version}")
    lines.append(f"")

    # Summary counts
    lines.append(f"## Summary")
    lines.append(f"")
    lines.append(f"| Category | Reference | Current | Difference |")
    lines.append(f"|----------|-----------|---------|------------|")

    for category in ('classes', 'enums', 'functions', 'constants'):
        ref_count = len(ref.get(category, {}))
        cur_count = len(cur.get(category, {}))
        diff = cur_count - ref_count
        diff_str = f"+{diff}" if diff > 0 else str(diff)
        lines.append(f"| {category.title()} | {ref_count} | {cur_count} | {diff_str} |")

    lines.append(f"")

    # --- Missing Classes ---
    ref_classes = set(ref.get('classes', {}).keys())
    cur_classes = set(cur.get('classes', {}).keys())
    missing_classes = sorted(ref_classes - cur_classes)
    added_classes = sorted(cur_classes - ref_classes)

    total_missing_methods = 0
    total_missing_props = 0
    classes_with_issues = {}

    # Check common classes for method/property differences
    for cls_name in sorted(ref_classes & cur_classes):
        issues = compare_classes(ref['classes'][cls_name], cur['classes'][cls_name], cls_name)
        if issues:
            classes_with_issues[cls_name] = issues
            for issue_type, issue_data in issues:
                if issue_type == 'missing_methods':
                    total_missing_methods += len(issue_data)
                elif issue_type == 'missing_properties':
                    total_missing_props += len(issue_data)

    # --- Missing Enums ---
    ref_enums = set(ref.get('enums', {}).keys())
    cur_enums = set(cur.get('enums', {}).keys())
    missing_enums = sorted(ref_enums - cur_enums)
    added_enums = sorted(cur_enums - ref_enums)

    enums_with_issues = {}
    for enum_name in sorted(ref_enums & cur_enums):
        issues = compare_enums(ref['enums'][enum_name], cur['enums'][enum_name], enum_name)
        if issues:
            enums_with_issues[enum_name] = issues

    # --- Missing Functions ---
    ref_funcs = set(ref.get('functions', {}).keys())
    cur_funcs = set(cur.get('functions', {}).keys())
    missing_funcs = sorted(ref_funcs - cur_funcs)
    added_funcs = sorted(cur_funcs - ref_funcs)

    # --- Overall Impact Summary ---
    lines.append(f"## Impact Summary")
    lines.append(f"")
    lines.append(f"| Issue | Count |")
    lines.append(f"|-------|-------|")
    lines.append(f"| Missing classes | {len(missing_classes)} |")
    lines.append(f"| Missing top-level enums | {len(missing_enums)} |")
    lines.append(f"| Missing top-level functions | {len(missing_funcs)} |")
    lines.append(f"| Classes with method differences | {len(classes_with_issues)} |")
    lines.append(f"| Total missing methods (in existing classes) | {total_missing_methods} |")
    lines.append(f"| Total missing properties (in existing classes) | {total_missing_props} |")
    lines.append(f"| Enums with value differences | {len(enums_with_issues)} |")
    lines.append(f"| Added classes (new in nanobind) | {len(added_classes)} |")
    lines.append(f"| Added top-level functions (new) | {len(added_funcs)} |")
    lines.append(f"")

    # --- Missing Classes Detail ---
    if missing_classes:
        lines.append(f"## Missing Classes ({len(missing_classes)})")
        lines.append(f"")
        lines.append(f"These classes exist in the Cython build but are missing from the nanobind build:")
        lines.append(f"")
        for cls in missing_classes:
            ref_cls = ref['classes'][cls]
            method_count = len(ref_cls.get('methods', {}))
            lines.append(f"- **{cls}** ({method_count} methods)")
        lines.append(f"")

    # --- Missing Top-Level Enums ---
    if missing_enums:
        lines.append(f"## Missing Top-Level Enums ({len(missing_enums)})")
        lines.append(f"")
        for e in missing_enums:
            ref_e = ref['enums'][e]
            val_count = len(ref_e.get('values', {}))
            lines.append(f"- **{e}** ({val_count} values)")
        lines.append(f"")

    # --- Missing Top-Level Functions ---
    if missing_funcs:
        lines.append(f"## Missing Top-Level Functions ({len(missing_funcs)})")
        lines.append(f"")
        for f in missing_funcs:
            sig = ref['functions'][f].get('signature', '')
            lines.append(f"- `{f}` {sig or ''}")
        lines.append(f"")

    # --- Class-Level Differences ---
    if classes_with_issues:
        lines.append(f"## Class-Level Differences ({len(classes_with_issues)} classes)")
        lines.append(f"")

        for cls_name, issues in sorted(classes_with_issues.items()):
            lines.append(f"### {cls_name}")
            lines.append(f"")
            for issue_type, issue_data in issues:
                if issue_type == 'missing_methods':
                    lines.append(f"**Missing methods:** {', '.join(f'`{m}`' for m in issue_data)}")
                elif issue_type == 'added_methods':
                    lines.append(f"**Added methods:** {', '.join(f'`{m}`' for m in issue_data)}")
                elif issue_type == 'missing_properties':
                    lines.append(f"**Missing properties:** {', '.join(f'`{p}`' for p in issue_data)}")
                elif issue_type == 'added_properties':
                    lines.append(f"**Added properties:** {', '.join(f'`{p}`' for p in issue_data)}")
                elif issue_type == 'missing_nested_enums':
                    lines.append(f"**Missing nested enums:** {', '.join(f'`{e}`' for e in issue_data)}")
                elif issue_type == 'missing_enum_values':
                    for ename, vals in issue_data.items():
                        lines.append(f"**Missing enum values in `{ename}`:** {', '.join(f'`{v}`' for v in vals)}")
                elif issue_type == 'static_changes':
                    lines.append(f"**Static/instance changes:** {', '.join(issue_data)}")
            lines.append(f"")

    # --- Enum-Level Differences ---
    if enums_with_issues:
        lines.append(f"## Enum-Level Differences ({len(enums_with_issues)} enums)")
        lines.append(f"")
        for enum_name, issues in sorted(enums_with_issues.items()):
            lines.append(f"### {enum_name}")
            for issue_type, issue_data in issues:
                if issue_type == 'missing_values':
                    lines.append(f"**Missing values:** {', '.join(f'`{v}`' for v in issue_data)}")
                elif issue_type == 'added_values':
                    lines.append(f"**Added values:** {', '.join(f'`{v}`' for v in issue_data)}")
                elif issue_type == 'value_changed':
                    lines.append(f"**Value changed:** {issue_data}")
            lines.append(f"")

    # --- Added Classes ---
    if added_classes:
        lines.append(f"## Added Classes ({len(added_classes)}) - New in Nanobind")
        lines.append(f"")
        lines.append(f"These classes exist in the nanobind build but NOT in the Cython build:")
        lines.append(f"")
        for cls in added_classes:
            cur_cls = cur['classes'][cls]
            method_count = len(cur_cls.get('methods', {}))
            lines.append(f"- **{cls}** ({method_count} methods)")
        lines.append(f"")

    # --- Added Functions ---
    if added_funcs:
        lines.append(f"## Added Top-Level Functions - New in Nanobind")
        lines.append(f"")
        for f in added_funcs:
            lines.append(f"- `{f}`")
        lines.append(f"")

    # --- Added Enums ---
    added_enum_list = sorted(cur_enums - ref_enums)
    if added_enum_list:
        lines.append(f"## Added Enums - New in Nanobind")
        lines.append(f"")
        for e in added_enum_list:
            lines.append(f"- `{e}`")
        lines.append(f"")

    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description='Compare pyopenms API JSON files')
    parser.add_argument('reference', help='Reference API JSON (e.g., Cython 3.5)')
    parser.add_argument('current', help='Current API JSON (e.g., nanobind 3.6)')
    parser.add_argument('--output', '-o', default=None, help='Output report file (default: stdout)')
    args = parser.parse_args()

    ref = load_api(args.reference)
    cur = load_api(args.current)

    report = generate_report(ref, cur)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(report)
        print(f"Report written to {args.output}", file=sys.stderr)
    else:
        print(report)


if __name__ == '__main__':
    main()

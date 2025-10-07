# How to Use the Issue Assessment

This guide explains how to use the issue assessment reports to triage OpenMS GitHub issues.

## Quick Start

1. **Read the Executive Summary**: Start with `ISSUE_ASSESSMENT_REPORT.md`
2. **Review Specific Issues**: Check `SPECIFIC_ISSUES_TO_REVIEW.md` for actionable items
3. **Use the Analysis Scripts**: Scripts in `/tmp` for future assessments

## Files Created

### Main Reports

#### `ISSUE_ASSESSMENT_REPORT.md`
- Executive summary of findings
- Statistical overview
- Pattern analysis
- Recommendations for issue triage

#### `SPECIFIC_ISSUES_TO_REVIEW.md`
- Detailed list of specific issues organized by priority
- Evidence for each assessment
- Recommended actions
- Direct links to GitHub issues

### Analysis Scripts (in /tmp)

#### `practical_issue_analyzer.py`
- Analyzes CHANGELOG for issue references
- Categorizes fixes by type
- Generates reports

**Usage:**
```bash
python3 /tmp/practical_issue_analyzer.py
# Outputs: /tmp/openms_issue_assessment.txt
#          /tmp/openms_issue_summary.json
```

#### `issue_assessment_script.py`
- Core assessment logic
- Can be imported and extended

#### `analyze_all_issues.py`
- Comprehensive analyzer when given issue data
- Generates markdown and JSON reports

**Usage:**
```bash
python3 /tmp/analyze_all_issues.py issues.json
```

## Workflow for Issue Triage

### Step 1: High Priority Issues (Likely Fixed)

1. Start with the "High Priority" section in `SPECIFIC_ISSUES_TO_REVIEW.md`
2. For each issue:
   - Read the evidence provided
   - Verify the fix exists (check CHANGELOG or test)
   - Add a comment: "This issue appears to be fixed in version X.Y. Verified by [evidence]. Closing."
   - Close the issue with appropriate label

Example issues to start with:
- #7271 (zlib-ng tests)
- #7112 (SageAdapter)
- #7171 (ARM64 test)

### Step 2: Medium Priority Issues (Possibly Fixed)

1. Review "Medium Priority" section
2. For each issue:
   - Request the original reporter to retest on current version
   - If no response in 30 days, close with: "Closing due to inactivity. Please reopen if still relevant in current version."
   - If confirmed fixed, close
   - If still broken, keep open and add details

### Step 3: Questions and Support

1. Identify question-type issues
2. Answer the question
3. Add: "Moving this to GitHub Discussions. Closing this issue."
4. Create equivalent discussion if needed
5. Close issue

### Step 4: Feature Requests

1. Review feature requests
2. Options:
   - If implemented: Close with reference to implementation
   - If planned: Keep open with "sprint candidate" label
   - If not planned: Close with explanation or convert to discussion
   - If needs design: Label as "needs design decision"

### Step 5: Code Review Needed

1. Issues marked "Requires Code Review"
2. Assign to appropriate maintainer
3. Schedule for code review
4. Fix or close based on findings

## Using the JSON Summary

The file `/tmp/openms_issue_summary.json` contains:

```json
{
  "fixed_issue_numbers": [23, 24, 27, ...],
  "total_fixed": 885,
  "changelog_categories": {
    "bugs": 411,
    "features": 766,
    ...
  }
}
```

**Use this to:**
- Batch process issue closures
- Track progress in spreadsheet
- Generate statistics
- Create GitHub API scripts

## Automation Example

You can create a script to check issues against the JSON:

```python
import json
import requests

# Load our analysis
with open('/tmp/openms_issue_summary.json') as f:
    data = json.load(f)

fixed_numbers = set(data['fixed_issue_numbers'])

# Check each open issue
for issue_num in fixed_numbers:
    # Use GitHub API to check if still open
    # If open, post comment and close
    pass
```

## Regular Maintenance

### Monthly
1. Re-run `practical_issue_analyzer.py` after releases
2. Review new issues against patterns
3. Update triage priorities

### Quarterly
1. Full audit of open issues
2. Update sprint candidates
3. Review feature request backlog
4. Archive very old issues

### After Major Releases
1. Request retesting for old version issues
2. Close issues explicitly fixed in release
3. Update documentation links

## Tips for Effective Triage

### DO:
- ✅ Be polite and appreciative of reporters
- ✅ Provide clear evidence when closing
- ✅ Link to CHANGELOG or commit when closing
- ✅ Thank users who verify fixes
- ✅ Convert questions to discussions
- ✅ Label issues consistently

### DON'T:
- ❌ Close issues without explanation
- ❌ Assume old issues are fixed without verification
- ❌ Leave users without response
- ❌ Close feature requests without explanation
- ❌ Ignore crash reports or data loss bugs

## Issue Labels to Use

When closing issues, add appropriate labels:

- `fixed`: Issue was fixed
- `duplicate`: Duplicate of another issue
- `wontfix`: Won't be fixed (explain why)
- `invalid`: Not actually an issue
- `question`: User question (move to discussions)
- `obsolete`: No longer relevant

## Priority Classification

| Priority | Criteria | Action Timeline |
|----------|----------|----------------|
| Critical | Crash, data loss, security | Immediate |
| High | Blocks users, incorrect results | Within week |
| Medium | Inconvenience, poor UX | Within month |
| Low | Enhancement, documentation | When time permits |

## Getting Help

If unsure about an issue:
1. Ask in OpenMS developer chat/forum
2. Mention relevant maintainer
3. Add "needs triage" label
4. Don't close - let others review

## Success Metrics

Track these metrics over time:
- Number of open issues (should decrease)
- Average issue age (should decrease)
- Issue response time (should be < 7 days)
- Percentage with labels (should be 100%)
- Number of stale issues (should be minimal)

## Example Comment Templates

### Closing Fixed Issue
```
This issue appears to be fixed in OpenMS X.Y.Z.

Evidence:
- CHANGELOG entry: [link]
- Commit: [link]
- Verified by testing [scenario]

Thank you for reporting this! Closing as fixed.
```

### Requesting Retest
```
This issue was reported for OpenMS X.X. We've made significant improvements since then.

Could you please retest with the current version (X.Y.Z) and let us know if the issue still exists?

If we don't hear back in 30 days, we'll close this issue. You can always reopen it or create a new one if the problem persists.
```

### Moving to Discussions
```
Thank you for your question! This is more suitable for GitHub Discussions rather than an issue.

I've answered your question there: [link]

Closing this issue and moving the conversation to discussions.
```

### Closing Feature Request
```
Thank you for the feature request!

[After review, we've decided to implement this in version X.Y / not implement because...]

[If implemented: Closing as completed]
[If not: Closing as wont-fix, but feel free to submit a PR if you'd like to contribute this feature]
```

## Questions?

If you have questions about this assessment or need help with triage:
1. Review the main reports first
2. Check the CHANGELOG for context
3. Ask in OpenMS developer channels
4. Create a meta-issue for triage process improvements

---

**Document Version**: 1.0  
**Last Updated**: 2025-01-08  
**Author**: GitHub Copilot Issue Assessment

---
description: Check OpenMS nightly build status, recent pull requests, and open issues. Gathers actionable items and prioritizes them.
allowed-tools: Bash(gh *), Bash(git *)
---

# OpenMS Status Dashboard

You are a project status analyst for the OpenMS project (https://github.com/OpenMS/OpenMS). Gather data from GitHub, analyze it, and present a prioritized action plan.

## Data Collection

Run ALL of the following `gh` commands in parallel to gather data:

### 1. Nightly Build Status
Check the most recent runs of key CI workflows on the `nightly` branch:

```bash
gh run list --repo OpenMS/OpenMS --branch nightly --limit 5 --json name,status,conclusion,createdAt,url,headBranch
```

Also check the main CI on `develop`:
```bash
gh run list --repo OpenMS/OpenMS --branch develop --workflow "openms-ci-full" --limit 3 --json name,status,conclusion,createdAt,url
```

If any workflow has `conclusion: "failure"`, drill into the failed jobs:
```bash
gh run view <run-id> --repo OpenMS/OpenMS --json jobs --jq '.jobs[] | select(.conclusion == "failure") | {name, conclusion, url}'
```

### 2. Recent Pull Requests
Fetch open PRs sorted by most recent activity:

```bash
gh pr list --repo OpenMS/OpenMS --state open --limit 20 --json number,title,author,createdAt,updatedAt,labels,reviewDecision,isDraft,url,reviewRequests,mergeable,headRefName
```

Also check recently merged PRs for context:
```bash
gh pr list --repo OpenMS/OpenMS --state merged --limit 5 --json number,title,mergedAt,author,url
```

### 3. Open Issues
Fetch open issues, prioritizing labeled ones:

```bash
gh issue list --repo OpenMS/OpenMS --state open --limit 30 --json number,title,labels,createdAt,updatedAt,author,url,comments --sort updated
```

Check for critical/blocker issues specifically:
```bash
gh issue list --repo OpenMS/OpenMS --state open --label critical,blocker,major --json number,title,labels,createdAt,url
```

### 4. New Issues
Fetch issues created in the last 14 days to surface newly reported bugs and feature requests:

```bash
gh issue list --repo OpenMS/OpenMS --state open --limit 20 --json number,title,labels,createdAt,author,url,comments --sort created
```

Filter the results to only those created within the last 14 days. These represent fresh reports that may need triage, labeling, or assignment.

## Analysis & Presentation

After collecting all data, present findings in this structure:

### Status Overview
A brief 2-3 line summary of overall project health.

### Nightly Build Status
- Show a table of recent nightly workflow runs with status indicators
- Use checkmarks/X marks for pass/fail
- If failures exist, list the specific failed jobs and platforms

### Action Items (Prioritized)

Organize actions into priority tiers:

**P0 - Critical (fix immediately)**
- Nightly build failures (broken builds block all development)
- Issues labeled `critical` or `blocker`
- PRs with failing CI that were previously approved

**P1 - High (address this week)**
- Issues labeled `major`
- PRs awaiting review for >7 days
- PRs with review changes requested but no recent activity
- Flaky or intermittent CI failures

**P2 - Medium (address soon)**
- PRs awaiting review for >3 days
- Open issues with recent activity/comments
- Draft PRs that have been open >14 days without updates

**P3 - Low (backlog)**
- Stale PRs/issues with no activity >30 days
- Enhancement requests without assignees

For each action item, include:
- A direct link to the PR/issue/workflow
- A one-line description of what needs to happen
- Who is involved (author, reviewer, assignee)

### New Issues (Last 14 Days)
- List newly created issues in a table with: issue number (linked), title, author, labels, and age
- Flag any that are unlabeled or unassigned — these need triage
- Highlight bug reports that may relate to recent PRs or nightly failures

### Recent Activity
- List the last 5 merged PRs as a quick changelog

### Recommended Actions
After all analysis, produce a numbered list of concrete recommended actions ordered by impact. Each action should:
- Start with a verb (e.g., "Fix", "Review", "Triage", "Merge", "Investigate")
- Reference the specific PR/issue/workflow by number or link
- Be specific enough to act on immediately (not vague like "look into CI")
- Include who should take the action if known

Example:
1. **Fix** nightly Linux build failure in #1234 — `FileHandler_test` segfault on Ubuntu 22.04
2. **Review** PR #5678 (open 9 days, no reviewer assigned) — adds MzTab export support
3. **Triage** new issue #9012 — user-reported crash in `FeatureFinderMetabo`, unlabeled

Aim for 5-10 actions. This list should serve as a ready-made to-do list for the maintainer.

## Output Format
- Use GitHub-flavored markdown tables where appropriate
- Keep the output scannable — use bold for key items, links for references
- End with the recommended actions list as the final section so it's easy to find

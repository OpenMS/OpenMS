#!/usr/bin/env python3
"""Sync the OpenDIA step board (BOARD.md) to GitHub Project OpenMS/38.

Source of truth: BOARD.md, sections "## R-WPn: ..." (the "## Now" table repeats rows and is not synced).
- Each R-WP issue (#10848..#10853) is added to the project, and its Status is set from WP_STATUS below.
- Each step row becomes one draft item titled "R-WPn: <step>". Its body has the BOARD.md state and the record link.
  An existing draft with the same title is updated, not duplicated. Drafts that are no longer in BOARD.md are left
  untouched and listed, or archived with --archive-stale (archiving is reversible in the project UI).
- State mapping to the project Status field:
    Done* -> Done, In progress / Partial -> In progress, Next -> Ready, Blocked / Open -> Backlog, Stopped -> Done.
Usage: board_sync.py [--dry-run] [--archive-stale]. Needs `gh` with the `project` scope.
"""
import json, re, subprocess, sys, pathlib

OWNER, NUMBER, REPO = 'OpenMS', '38', 'OpenMS/OpenMS'
HERE = pathlib.Path(__file__).resolve().parent
BOARD = HERE / 'BOARD.md'
DOC_URL = 'https://github.com/OpenMS/OpenMS/blob/docs/opendia-wp00-design-plan/doc/design/opendia-laptop-first/'
WP_ISSUES = {'R-WP1': 10848, 'R-WP2': 10849, 'R-WP3': 10850, 'R-WP4': 10851, 'R-WP5': 10852, 'R-WP6': 10853}
WP_STATUS = {'R-WP1': 'In review', 'R-WP2': 'In review', 'R-WP3': 'In progress', 'R-WP4': 'In progress',
             'R-WP5': 'In progress', 'R-WP6': 'Backlog'}
DRY = '--dry-run' in sys.argv
ARCHIVE = '--archive-stale' in sys.argv


def gh(*args):
    out = subprocess.run(['gh', *args], check=True, capture_output=True, text=True, timeout=120).stdout
    return json.loads(out) if out.strip().startswith(('{', '[')) else out


def status_of(state):
    s = re.sub(r'[*]', '', state).strip().lower()
    if s.startswith('done') or s.startswith('stopped'):
        return 'Done'
    if s.startswith('in progress') or s.startswith('partial'):
        return 'In progress'
    if s.startswith('next'):
        return 'Ready'
    return 'Backlog'


def parse_board():
    rows, wp = [], None
    for line in BOARD.read_text().splitlines():
        m = re.match(r'## (R-WP\d)\b', line)
        if m:
            wp = m.group(1)
            continue
        if line.startswith('## '):
            wp = None
        if not wp or not line.startswith('|') or line.startswith('|---') or line.startswith('| Step'):
            continue
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        if len(cells) < 3:
            continue
        step, state, record = cells[0], cells[1], cells[2]
        # Relative markdown links in the record cell point into the docs tree.
        record = re.sub(r'\]\((?!http)([^)]+)\)', lambda m: '](' + DOC_URL + m.group(1) + ')', record)
        title = f'{wp}: ' + re.sub(r'`', '', re.sub(r'\*\*', '', step))
        if len(title) > 250:
            title = title[:247] + '...'
        body = (f'**State (BOARD.md):** {state}\n\n**Record:** {record or "none yet"}\n\n'
                f'**Work package:** #{WP_ISSUES[wp]}\n\nSynced from [BOARD.md]({DOC_URL}BOARD.md) by `board_sync.py`.')
        rows.append({'wp': wp, 'title': title, 'body': body, 'status': status_of(state)})
    return rows


def main():
    project_id = gh('project', 'view', NUMBER, '--owner', OWNER, '--format', 'json')['id']
    fields = gh('project', 'field-list', NUMBER, '--owner', OWNER, '--format', 'json')['fields']
    status_field = next(f for f in fields if f['name'] == 'Status')
    option = {o['name']: o['id'] for o in status_field['options']}
    items = gh('project', 'item-list', NUMBER, '--owner', OWNER, '--limit', '500', '--format', 'json')['items']
    by_issue = {i['content'].get('number'): i for i in items if i['content'].get('type') == 'Issue'}
    drafts = {i['title']: i for i in items if i['content'].get('type') == 'DraftIssue'}
    # item-list omits the DraftIssue content id (DI_...), which item-edit needs for title/body.
    draft_ids, cursor = {}, None
    while True:
        after = f', after: "{cursor}"' if cursor else ''
        q = ('{organization(login: "%s") {projectV2(number: %s) {items(first: 100%s) {pageInfo {hasNextPage endCursor} '
             'nodes {id content {... on DraftIssue {id}}}}}}}' % (OWNER, NUMBER, after))
        page = gh('api', 'graphql', '-f', f'query={q}')['data']['organization']['projectV2']['items']
        draft_ids.update({n['id']: n['content']['id'] for n in page['nodes'] if n.get('content') and n['content'].get('id')})
        if not page['pageInfo']['hasNextPage']:
            break
        cursor = page['pageInfo']['endCursor']

    def set_status(item_id, name):
        if DRY:
            return
        gh('project', 'item-edit', '--id', item_id, '--project-id', project_id, '--field-id', status_field['id'],
           '--single-select-option-id', option[name])

    log = []
    for wp, num in WP_ISSUES.items():
        item = by_issue.get(num)
        if item is None:
            if not DRY:
                item = gh('project', 'item-add', NUMBER, '--owner', OWNER, '--url',
                          f'https://github.com/{REPO}/issues/{num}', '--format', 'json')
            log.append(f'add issue #{num} ({wp})')
        if item is not None and item.get('status') != WP_STATUS[wp]:
            set_status(item['id'], WP_STATUS[wp])
            log.append(f'status #{num} -> {WP_STATUS[wp]}')
    seen = set()
    for r in parse_board():
        seen.add(r['title'])
        d = drafts.get(r['title'])
        if d is None:
            if not DRY:
                d = gh('project', 'item-create', NUMBER, '--owner', OWNER, '--title', r['title'], '--body', r['body'],
                       '--format', 'json')
            log.append(f'create draft "{r["title"]}" [{r["status"]}]')
        else:
            if not DRY and d['content'].get('body') != r['body']:
                gh('project', 'item-edit', '--id', draft_ids[d['id']], '--title', r['title'], '--body', r['body'])
                log.append(f'update draft "{r["title"]}"')
        if d is not None and d.get('status') != r['status']:
            set_status(d['id'], r['status'])
            log.append(f'  status -> {r["status"]}')
    stale = [t for t in drafts if t.startswith('R-WP') and t not in seen]
    if stale and ARCHIVE:
        for t in stale:
            if not DRY:
                gh('project', 'item-archive', NUMBER, '--owner', OWNER, '--id', drafts[t]['id'])
            log.append(f'archive draft "{t}"')
        stale = []
    print('\n'.join(log))
    if stale:
        print('drafts not in BOARD.md (left untouched):\n  ' + '\n  '.join(stale))
    print(f'{"dry run: " if DRY else ""}{len(log)} actions')


if __name__ == '__main__':
    main()

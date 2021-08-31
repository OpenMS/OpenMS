from github_methods import *
from github import Github
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Create or update GitHub issue with identified unknown words.'
                                                 'Only correctly callable from GitHub-push action.')
    parser.add_argument('-t', '--token', type=str)
    parser.add_argument('-repo', '--repository', type=str)
    parser.add_argument('-br', '--branch', type=str)
    parser.add_argument('-c', '--commit', type=str)
    parser.add_argument('-f', '--full', action="store_true")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    g = Github(args.token)

    repo = g.get_repo(args.repository)

    branch = args.branch.split('/')[-1]
    title = f'Spell-Check Results - {args.repository.split("/")[0]}/{branch}'

    # Find out if issue already exists
    issue = [issue for issue in repo.get_issues() if issue.title == title]

    edited_files = False
    if not args.full:
        commit = repo.get_commit(args.commit)
        edited_files = {Path(file.filename) for file in commit.files}

    unknown_words = get_words(edited_files)

    if len(unknown_words) > 0:

        if len(issue) > 0:
            issue = issue[0]

            if not args.full:
                old_unknown_words = body_to_words(issue.body)
                for word, properties in unknown_words.items():
                    for file in list(old_unknown_words[word]['files'].keys()):
                        if file in edited_files and file not in unknown_words[word]['files']:
                            old_unknown_words[word]['files'].pop(file)
                    for file, lines in properties['files'].items():
                        old_unknown_words[word]['files'][file] = lines
                unknown_words = old_unknown_words

            issue_body = words_to_body(title, unknown_words)
            issue.edit(body=issue_body)
        else:
            issue_body = words_to_body(title, unknown_words)
            issue = repo.create_issue(title, issue_body)


if __name__ == '__main__':
    main()

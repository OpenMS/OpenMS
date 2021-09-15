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
    title = f'Spellcheck Results - {args.repository.split("/")[0]}/{branch}'

    # Find out if issue already exists
    issue = [issue for issue in repo.get_issues() if issue.title == title]

    # Default: run search only on edited files
    edited_files = False
    if not args.full:
        commit = repo.get_commit(args.commit)
        edited_files = {Path(file.filename) for file in commit.files}

    print('Searching words..')
    unknown_words = get_words(files_filter=edited_files)
    print(f'Search finished, {len(unknown_words)} unknown words were found!')
    update_vocab(repo, branch)

    if len(unknown_words) > 0:

        # Issue exists already
        if len(issue) > 0:
            issue = issue[0]

            if not args.full:
                print('Processing current GitHub Issue..')
                old_unknown_words = comments_to_words(issue.get_comments())
                print(f'GitHub Issue processed, {len(unknown_words)} words retained!')

                # Word got deleted in edited files
                print('Removing deleted words and occurences...')
                len_before_removing = len(unknown_words)
                removed_files = 0
                for word in list(old_unknown_words.keys()):
                    for file in list(old_unknown_words[word]['files'].keys()):
                        if file in edited_files:
                            if word in unknown_words:
                                if file not in unknown_words[word]['files']:
                                    old_unknown_words[word]['files'].pop(file)
                                    removed_files += 1
                                if len(old_unknown_words[word]['files']) == 0:
                                    old_unknown_words.pop(word)
                print(f'{len(old_unknown_words) - len_before_removing} removed words and '
                      f'{removed_files} files were removed!')

                # Updated, added words in edited files
                print('Transferring old words and properties...')
                for word, properties in unknown_words.items():
                    for file, lines in properties['files'].items():
                        old_unknown_words[word]['files'][file] = lines
                unknown_words = {key: old_unknown_words[key] for key in
                                 sorted(old_unknown_words.keys(), key=str.casefold)}
                print('Done!')
        else:
            # Create new issue
            issue = repo.create_issue(title, ' ', labels=['spellcheck'])
        print('Convert unknown words to commments...')
        comments = words_to_comments(unknown_words, repo, branch)
        print('Done!')

        print('Update issue...')
        update_issue(issue, title, comments, len(unknown_words))
    print('All tasks executed successfully')


if __name__ == '__main__':
    main()

from github_methods import *
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Process a Spell-Check Results - GitHub issue.'
                                                 'Only correctly callable from a GitHub-issue edit action.')
    parser.add_argument('-to', '--token', type=str)
    parser.add_argument('-repo', '--repository', type=str)
    parser.add_argument('-ti', '--title', type=str)
    parser.add_argument('-in', '--issue_number', type=int)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    g = Github(args.token)

    repo = g.get_repo(args.repository)
    issue = repo.get_issue(args.issue_number)
    branch = args.title.split('/')[-1]
    title = f'Spell-Check Results - {args.repository.split("/")[0]}/{branch}'

    unknown_words = comments_to_words(issue.get_comments())
    process_actions_github(unknown_words, repo, branch)

    new_vocabulary = list(vocabulary)
    new_vocabulary.sort()
    new_vocabulary_file = json.dumps(new_vocabulary, indent=4)

    vocabulary_file = repo.get_contents(PATH_VOCAB, branch)
    repo.update_file(PATH_VOCAB, 'Update words', new_vocabulary_file, vocabulary_file.sha, branch)

    # Update issue comments
    comments = words_to_comments(unknown_words)
    update_issue(issue, title, comments, len(unknown_words))


if __name__ == '__main__':
    main()

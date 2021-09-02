from github_methods import *
from github import Github
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

    reference = set_ref()

    unknown_words = comments_to_words(issue.get_comments())
    process_actions_github(reference, unknown_words, repo, branch)

    new_vocabulary_file_content = json.dumps(VOCABULARY, indent=4)

    vocabulary_file = repo.get_contents(VOCAB_PATH, branch)
    repo.update_file(VOCAB_PATH, 'Update words', new_vocabulary_file_content, vocabulary_file.sha, branch)

    # Update issue comments
    comments = words_to_comments(unknown_words)
    for comment in issue.get_comments():
        comment.delete()
    for comment in comments:
        issue.create_comment(comment)


if __name__ == '__main__':
    main()

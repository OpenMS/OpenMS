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

    vocabulary_file = repo.get_contents(PATH_VOCABULARY)
    vocabulary = json.loads(vocabulary_file.decoded_content)
    reference = set_ref(vocabulary)

    unknown_words = body_to_words(issue.body)
    process_actions_github(reference, unknown_words, repo, branch)

    new_vocabulary_file_content = json.dumps(vocabulary, indent=4)
    repo.update_file(PATH_VOCABULARY, 'Update words', new_vocabulary_file_content, vocabulary_file.sha, branch)

    # Update issue
    title = f'Spell-Check Results - {args.repository.split("/")[0]}/{branch}'
    issue_body = words_to_body(title, vocabulary, unknown_words)
    issue.edit(body=issue_body)


if __name__ == '__main__':
    main()

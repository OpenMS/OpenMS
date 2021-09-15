from utils import *
from github.Repository import Repository
from github import ContentFile, PaginatedList, Issue


INFORMATION = """
Please read the provided README.md in `tools/spellcheck` carefully before continuing.
Use one of the options for each word. 

Important: Check the "Confirm" comments first. The words are accepted by the spellcheck, 
but are not placed in the vocabulary yet.

Replacing the word will automatically put the replacement in the vocabulary.
Replace the _-underscore and check the box in the last replacing option to use a custom word.
When you are done, create another comment containing only the word "process".
This will trigger the processing and resolve any potential issues.
"""


def words_to_comments(unknown_words: Union[dict, defaultdict], repo: Repository, branch: str) -> list:
    """
    Convert unknown words to github issue body

    :param unknown_words: All unknown, unprocessed words
    :param repo: The GitHub repository
    :param branch: Current branch of repository
    :return: Issue body
    """

    comments = ['']
    for i_word, word in enumerate(sorted(unknown_words.keys(), key=str.casefold)):
        error, action, blob, files = unknown_words[word].values()
        word_block = f'[{i_word+1}] "**{word}**" in file(s):\n'
        if len(comments[-1]) >= 50000:
            comments.append('')
        for i_path, (path, lines) in enumerate(list(files.items())):
            if i_path >= 5:
                word_block += '`...`\n'
                break
            word_block += f'<details>\n<summary>{path}</summary>\n\n'
            content_file = repo.get_contents(str(path), branch)
            file = open(path)
            file_lines = file.readlines()
            for i_line, line in enumerate(str(lines)[1:-1].split(', ')):
                if i_line >= 5:
                    word_block += '`...`\n'
                    break
                word_block += f'({line}) [`{file_lines[int(line)-1][:-1]}`]({content_file.html_url}#L{line})\n'
            word_block += '</details>\n'

        if error != '':
            word_block += f'{error}\n'

        word_block += f'\n- [ ] Add "{word}" to vocabulary *({blob[0][1]}%)*'
        if len(blob) > 1:
            for blob_word, certainty in blob[1:]:
                word_block += f'\n- [ ] Replace "{word}" with: "{blob_word}" *({certainty}%)*'
        word_block += f'\n- [ ] Replace "{word}" with "_" *(custom)*'
        comments[-1] += word_block + '\n\n\n'
    return comments


def comments_to_words(comments: PaginatedList) -> defaultdict:
    """
    Convert github issue body to unknown words dictionary

    :param comments: All comments of the issue
    :return: Unknown words
    """
    unknown_words = defaultdict(lambda: {'error': '', 'action': {'replacement': '', 'vocabulary': ''},
                                         'blob': [], 'files': defaultdict(list)})

    for comment in comments:
        lines = comment.body.split('\n')
        i = 0
        while i < len(lines):
            # Word block begin
            if lines[i].startswith('['):
                word = lines[i].split('**')[1]

                path = ''
                while not lines[i].startswith('-'):

                    if lines[i].startswith('<summary>'):
                        path = lines[i][9:-10]

                    if lines[i].startswith('(') and path != '':
                        unknown_words[word]['files'][path].append(int(lines[i][1:].split(')')[0]))
                    i += 1

                # Vocabulary
                vocab_action = word if lines[i][3].casefold() == 'x' else ''
                bwords = [[lines[i].split('"')[1], float(lines[i].split('(')[-1].split(')')[0][:-1])]]

                i += 1

                # Replacing
                replace_action = ''
                while lines[i].startswith('-'):
                    choice = lines[i][3].casefold()
                    replacement = lines[i].split('"')[-2]
                    certainty = lines[i].split('(')[-1].split(')')[0]
                    if certainty != 'custom':
                        bwords.append([replacement, float(certainty[:-1])])
                    if choice == 'x':
                        if replace_action != '':
                            unknown_words[word]['error'] = f'Error: You cannot assign multiple actions to a word.'
                            replace_action = ''
                            break
                        else:
                            replace_action = replacement
                    i += 1

                unknown_words[word]['blob'] = bwords

                if vocab_action != '':
                    if replace_action != '':
                        unknown_words[word]['error'] = f'Error: You cannot assign multiple actions to a word.'
                    else:
                        unknown_words[word]['action']['vocabulary'] = vocab_action
                elif replace_action != '':
                    unknown_words[word]['action']['replacement'] = replace_action
                    unknown_words[word]['action']['vocabulary'] = replace_action
            i += 1
    return unknown_words


def update_issue(issue: Issue, title: str, comments: list, len_unknown_words: int):
    """
    Update the GitHub issue including replacing all comments

    :param issue: GitHub issue
    :param title: GitHub issue title
    :param comments: New comments
    :param len_unknown_words: Number of unknown words
    """
    body = f"---\n{title}\n---\n\n" \
           f"{INFORMATION}\n\n" \
           f"Found {len_unknown_words} unknown words!"
    issue.edit(body=body)

    issue.unlock()
    for comment in issue.get_comments():
        comment.delete()

    for n_comment, comment in enumerate(comments):
        issue.create_comment(comment)
    issue.lock('resolved')


def process_actions_github(unknown_words: Union[dict, defaultdict], repo: Repository, branch: str):
    """
    Process actions of unknown words in github

    :param unknown_words: All unknown, unprocessed words
    :param repo: GitHub repository
    :param branch: Current branch of repository
    """
    words = list(unknown_words.keys())
    edited_files = defaultdict(lambda: {'content': [], 'content_file': ContentFile})
    for word in words:
        error, action, blob, files = unknown_words[word].values()
        replacement, vocab_word = action.values()

        if replacement != '':
            vocabulary.add(replacement)

            # Build files from all replacements
            for path, lines in files.items():
                if path not in edited_files:
                    edited_files[path]['content_file'] = repo.get_contents(str(path), branch)
                    edited_files[path]['content'] = edited_files[path]['content_file'].decoded_content.decode()
                content = edited_files[path]['content']
                new_content = ''
                for n_line, line in enumerate(content.split('\n')):
                    if n_line+1 in lines:
                        line = line.replace(word, replacement)
                    new_content += f'{line}\n'
                edited_files[path]['content'] = new_content.strip() + '\n'

            unknown_words.pop(word)

        elif vocab_word != '':
            vocabulary.add(vocab_word)
            unknown_words.pop(vocab_word)

    # Process all files affected
    for path, properties in edited_files.items():
        repo.update_file(str(path), f'Replaced words', properties['content'].encode(),
                         properties['content_file'].sha, branch)


def update_vocab(repo, branch):
    """
    Update the GitHub vocabulary
    """
    new_vocabulary = list(vocabulary)
    new_vocabulary.sort()
    new_vocabulary_file = json.dumps(new_vocabulary, indent=4)

    vocabulary_file = repo.get_contents('.github/spellcheck/vocabulary.json', branch)
    repo.update_file('.github/spellcheck/vocabulary.json', 'Update words', new_vocabulary_file, vocabulary_file.sha,
                     branch)

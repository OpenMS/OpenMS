from spellcheck import *
from github.Repository import Repository
from github import ContentFile


INFORMATION = """
Please read the provided README.md in `tools/spellcheck` carefully before continuing.
State the replacement and or the vocabulary index by replacing the whitespace in the respective ` ` code-box.
Word replacements, that are assigned a vocabulary index, will be ignored, if the replacement already exists in
the vocabulary.
"""


def words_to_body(title: str, vocabulary: dict, unknown_words: defaultdict) -> str:
    """
    Convert unknown words to github issue body

    :param title: Title of issue body
    :param vocabulary: Vocabulary of all whitelisted words
    :param unknown_words: All unknown, unprocessed words
    :return: Issue body
    """
    issue_body = f'---\n{title}\n---\n'
    issue_body += f'{INFORMATION}\n\n'
    issue_body += get_vocab_keys(vocabulary, '::', '......', '>')
    issue_body += '\n\n'

    for i, (word, properties) in enumerate(unknown_words.items()):
        issue_body += f'[{i + 1}] "**{word}**" in file(s):\n'
        for file, lines in list(properties['files'].items()):
            issue_body += f'`{file}: {str(lines)[1:-1]}`\n'
        if properties['error'] != '':
            issue_body += f'{properties["error"]}\n'
        issue_body += f'\n- Replace "{word}" with: ` `\n'
        issue_body += '- Add to Vocabulary: ` `\n\n\n'

    return issue_body


def body_to_words(issue_body: str) -> [defaultdict, str]:
    """
    Convert github issue body to unknown words dictionary

    :param issue_body: Issue body
    :return: Unknown words
    """
    unknown_words = defaultdict(lambda: {'error': '', 'action': {'replacement': '', 'vocab_index': ''},
                                         'files': defaultdict(list)})
    lines = issue_body.split('\n')
    i = 0
    while i < len(lines):
        if lines[i].startswith('['):
            word = lines[i].split('**')[1]
            i += 1
            while lines[i].startswith('`'):
                file, file_lines = lines[i][1:-1].split(':')
                unknown_words[word]['files'][Path(file)] = [int(l) for l in file_lines.strip(' ').strip('`').split(',')]
                i += 1
            while not lines[i].startswith('-'):
                i += 1
            replace_action = lines[i].split('`')[-2].strip()
            vocab_action = lines[i + 1].split('`')[-2].strip()
            if vocab_action != '':
                try:
                    vocab_action = int(vocab_action)
                    unknown_words[word]['action']['vocab_index'] = vocab_action
                except ValueError:
                    unknown_words[word]['error'] = f'Input "{vocab_action}" for "{word}" is invalid'
            if replace_action != '':
                unknown_words[word]['action']['replacement'] = replace_action
        i += 1
    return unknown_words


def process_actions_github(reference: dict, unknown_words: defaultdict, repo: Repository, branch: str):
    """
    Process actions of unknown words in github

    :param reference: Reference dictionary for each list containing header of vocabulary
    :param unknown_words: All unknown, unprocessed words
    :param repo: Github repository
    :param branch: Current branch of repository
    """
    words = list(unknown_words.keys())
    files = defaultdict(lambda: {'content': [], 'content_file': ContentFile})
    for word in words:
        properties = unknown_words[word]
        replacement, vocab_index = properties['action'].values()
        in_vocab = any([replacement in set(words) for words in reference.values()])

        if replacement != '':
            if vocab_index == '':
                unknown_words[replacement] = unknown_words[word]

            for path, lines in properties['files'].items():
                if path not in files:
                    files[path]['content_file'] = repo.get_contents(str(path))
                    files[path]['content'] = files[path]['content_file'].decoded_content.decode()
                content = files[path]['content']
                new_content = ''
                for n_line, line in enumerate(content.split('\n')):
                    if n_line in lines:
                        line = line.replace(word, replacement)
                    new_content += f'{line}\n'
                files[path]['content'] = new_content.strip()

        if vocab_index != '':
            if vocab_index < 1 or vocab_index > len(reference):
                unknown_words[word]['error'] = f'Invalid index, number must be 1-{len(reference)}.'
                unknown_words[word]['action']['vocab_index'] = ''
            if not in_vocab:
                reference[vocab_index - 1].append(word if replacement == '' else replacement)
                reference[vocab_index - 1].sort(key=str.casefold)

        if replacement != '' or vocab_index != '':
            unknown_words.pop(word)

    # Process all files affected
    for path, properties in files.items():
        repo.update_file(str(path), f'Replaced words', properties['content'].encode(),
                         properties['content_file'].sha, branch)

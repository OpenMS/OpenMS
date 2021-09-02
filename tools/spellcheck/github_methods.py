from spellcheck import *
from github.Repository import Repository
from github import ContentFile, PaginatedList


def words_to_comments(unknown_words: defaultdict) -> list:
    """
    Convert unknown words to github issue body

    :param unknown_words: All unknown, unprocessed words
    :return: Issue body
    """
    comments = []
    comment = ''
    for i, word in enumerate(sorted(unknown_words.keys(), key=str.casefold)):
        properties = unknown_words[word]
        if len(comment) > 50000:
            comments.append(str(comment))
            comment = ''
        comment += f'[{i + 1}] "**{word}**" in file(s):\n'
        for file, lines in list(properties['files'].items()):
            comment += f'`{file}: {str(lines)[1:-1]}`\n'
        if properties['error'] != '':
            comment += f'{properties["error"]}\n'
        comment += f'\n- Replace "{word}" with: ` `\n'
        comment += '- Add to Vocabulary: ` `\n\n\n'
    comments.append(comment)
    return comments


def comments_to_words(comments: PaginatedList) -> defaultdict:
    """
    Convert github issue body to unknown words dictionary

    :param comments: All comments of the issue
    :return: Unknown words
    """
    unknown_words = defaultdict(lambda: {'error': '', 'action': {'replacement': '', 'vocab_index': ''},
                                         'files': defaultdict(list)})

    for comment in comments:
        lines = comment.body.split('\n')
        i = 0
        while i < len(lines):
            if lines[i].startswith('['):
                word = lines[i].split('**')[1]
                i += 1
                while lines[i].startswith('`'):
                    file, file_lines = lines[i][1:-1].split(':')
                    unknown_words[word]['files'][Path(file)] = [int(l) for l in
                                                                file_lines.strip(' ').strip('`').split(',')]
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

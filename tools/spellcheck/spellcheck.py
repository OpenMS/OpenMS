from pathlib import Path
import os
from collections import defaultdict
import json
import re
from typing import Union

PATH_VOCABULARY = 'vocabulary.json'
PATH_UNKNOWN_WORDS = 'unknown_words.json'
PATH_RULES = 'rules.json'


def clear():
    """Clears stdout"""
    os.system('clear')


def load_json(path: str):
    """
    Load json file from str

    :param path: str of file path
    """
    with open(path, 'r') as file:
        return json.load(file)


def write_json(obj: dict, path: str):
    """
    Write dict to .json file

    :param obj: Dictionary to save
    :param path: str of file path
    """
    with open(path, 'w') as file:
        json.dump(obj, file, indent=4)


def build_file_list(rules: dict) -> set:
    """
    Build list of all files to search, based on the rules provided.

    :param rules: Dictionary to save
    :return: set of included files
    """
    incl_ext = set(rules['include']['extensions'])
    incl = set()
    excl_ext = set(rules['exclude']['extensions'])
    for included in rules['include']['locations']:
        for path in Path(included).rglob("*"):
            if path.is_file():
                file_ext = path.suffix
                if (len(incl_ext) > 0 and file_ext in incl_ext) or len(incl_ext) == 0:
                    if file_ext not in excl_ext \
                            and not any([path.is_relative_to(excl_loc) for excl_loc in rules['exclude']['locations']]):
                        incl.add(path)
    return incl


def flatten_vocab(vocabulary) -> set:
    """
    Flatten a nested dictionary.

    :return: Flattened dictionary
    """
    set_out = set()

    def _inner(section):
        for key, values in section.items():
            set_out.add(key)
            if type(values) == dict:
                _inner(values)
            else:
                for entry in values:
                    set_out.add(entry)

    _inner(vocabulary)
    return set_out


def set_ref(vocabulary: dict) -> dict:
    """
    Reference the nested vocabulary correct_words

    :param vocabulary: Vocabulary of all whitelisted words
    :return: Reference dictionary for each list containing header of vocabulary
    """
    i = 0
    ref = dict()

    def _inner(dct):
        nonlocal i
        for key, value in dct.items():
            if type(value) == dict:
                _inner(value)
            else:
                ref[i] = value
                i += 1

    _inner(vocabulary)
    return ref


def get_vocab_keys(vocabulary, header: str = '', indent: str = ' ', prefix: str = '') -> str:
    """
    Return a pretty print string from all keys from vocabulary.

    :param vocabulary: Vocabulary of all whitelisted words
    :param header: Header string prefix
    :param indent: Indentation string prefix
    :param prefix: Prefix used for all headers equally
    :return: str representation of all vocabulary headers
    """
    i = 1
    printable = []

    def _inner(dct: dict, indt: str):
        nonlocal i
        for key, value in dct.items():
            if type(value) == dict:
                printable.append(f'{prefix}{indt}{header * (i // 10 + 1)} {key}\n')
                _inner(value, indt + indent)
            else:
                printable.append(f'{prefix}{indt} {i} {key}\n')
                i += 1

    _inner(vocabulary, '')
    return ''.join(printable)


COMMENT_TYPES = load_json('comment_types.json')


def get_words(vocabulary: dict, files_filter: Union[set, bool] = False, verbose: bool = False):
    """
    Find all valid words from all files defined by rules.json.

    :param vocabulary: Vocabulary of all whitelisted words
    :param files_filter: Filter for included_files
    :param verbose: Verbosity of word search in files
    :return: All valid words from all included files
    """
    unknown_words = defaultdict(lambda: {'error': '', 'action': {'replacement': '', 'vocab_index': ''},
                                         'files': defaultdict(list)})
    flat_vocab = flatten_vocab(vocabulary)
    rules = load_json(PATH_RULES)
    included_files = build_file_list(rules)
    errors = []
    file_list = files_filter.intersection(included_files) if files_filter else included_files
    pattern = rules['include']['pattern']

    def _search_file():
        try:
            file = open(path, 'r')
            txt = file.read()
            file_ext = path.suffix
            block_comment = ''
            for n_line, line in enumerate(txt.split('\n')):
                line_comment = ''
                for txt_block in line.split():
                    if file_ext in COMMENT_TYPES:
                        if txt_block == block_comment:
                            block_comment = ''
                        elif txt_block in {ext.split()[0] for ext in COMMENT_TYPES[file_ext]['block']}:
                            block_comment = {(e := ext.split())[0]: e[1] for ext in
                                             COMMENT_TYPES[file_ext]['block']}[txt_block]
                            continue
                        elif txt_block in COMMENT_TYPES[file_ext]['line']:
                            line_comment = txt_block
                            continue
                    if block_comment != '' or line_comment != '' or file_ext not in COMMENT_TYPES:
                        for word in re.findall(pattern, txt_block):
                            if word not in flat_vocab and \
                                    all([re.fullmatch(p, word) is None for p in rules['exclude']['patterns']]):
                                unknown_words[word]['files'][str(os.path.relpath(path, 'Spellcheck'))[3:]].append(n_line)
            file.close()
        except UnicodeDecodeError:
            errors.append(f'{UnicodeDecodeError}: {path}')

    if verbose:
        from progress.bar import Bar
        with Bar('Searching files..', max=len(file_list)) as bar:
            for path in file_list:
                _search_file()
                bar.next()
    else:
        for path in file_list:
            _search_file()
    return unknown_words

from pathlib import Path
import os
from collections import defaultdict
import json
import re
from typing import Union
from progress.bar import Bar, IncrementalBar
from textblob import Word


SP_DIR = os.path.dirname(os.path.realpath(__file__)) + '/'
PATH_VOCABULARY = Path(SP_DIR + 'vocabulary.json')
PATH_UNKNOWN_WORDS = Path(SP_DIR + 'unknown_words.json')


def clear():
    """Clears stdout"""
    os.system('clear')


def load_json(path: Path):
    """
    Load json file from str

    :param path: str of file path
    """
    with open(path, 'r') as file:
        return json.load(file)


def create_patterns():
    COMMENT_TYPES = load_json(Path(SP_DIR + 'comment_types.json'))
    comment_patterns = {}
    for ext, properties in COMMENT_TYPES.items():
        neg_lookbehind = '(?<![\'"])'
        neg_lookahead = '(?![\'"])'

        neg_lc_lbehind = ''.join([f'(?<!{re.escape(lc)})' for lc in properties['line']])
        line_comments = '(' + '|'.join([f'{re.escape(lc)}' for lc in properties['line']]) + ')'
        pattern = neg_lookbehind + line_comments + '.*' + neg_lookahead + '|'

        for start, end in properties['block']:
            pattern += f'{neg_lookbehind}{neg_lc_lbehind}{re.escape(start)}(.|\n)*?{re.escape(end)}{neg_lookahead}|'
        comment_patterns[ext] = pattern[:-1]
    return comment_patterns


RULES = load_json(Path(SP_DIR + 'rules.json'))
COMMENT_PATTERNS = create_patterns()
vocabulary = set(load_json(Path(SP_DIR + 'vocabulary.json')))


def save_json(obj: Union[dict, list], path: Path):
    """
    Write dict or list to .json file

    :param obj: Dictionary to save
    :param path: str of file path
    """
    with open(path, 'w') as file:
        json.dump(obj, file, indent=4)


def save_vocab():
    """
    Write vocabulary to vocabulary.json file
    """
    new_vocab = list(vocabulary)
    new_vocab.sort()
    save_json(new_vocab, PATH_VOCABULARY)


def build_file_list(rel_path: str = '', files_filter: Union[set, bool] = False) -> set:
    """
    Build list of all files to search, based on the rules provided.

    :param rel_path: Relative path prefix (needed for local execution)
    :param files_filter: Filter for included_files
    :return: set of included files
    """
    incl_ext = set(RULES['include']['extensions'])
    incl = set()
    excl_ext = set(RULES['exclude']['extensions'])
    for included in RULES['include']['locations']:
        for path in Path(rel_path + included).rglob("*"):
            if path.is_file():
                file_ext = path.suffix
                if (len(incl_ext) > 0 and file_ext in incl_ext) or len(incl_ext) == 0:
                    if file_ext not in excl_ext \
                            and not any([path.is_relative_to(Path(rel_path + excl_loc))
                                         for excl_loc in RULES['exclude']['locations']]):
                        incl.add(path)
    file_list = files_filter.intersection(incl) if files_filter else incl
    return file_list


def get_words(rel_path: str = '', files_filter: Union[set, bool] = False) -> defaultdict:
    """
    Find all valid words from all files defined by rules.json.

    :param rel_path: Relative path prefix (needed for local execution)
    :param files_filter: Filter for included_files
    :return: All valid words from all included files
    """
    unknown_words = defaultdict(lambda: {'error': '', 'action': {'replacement': '', 'vocabulary': ''},
                                         'blob': [], 'files': defaultdict(list)})
    file_list = build_file_list(rel_path, files_filter)
    include_pattern = RULES['include']['pattern']
    errors = []

    def _search_file():
        """
        Search all words within a file.
        """

        def _search_block(text_block):
            """
            Search all words within a text block.
            """
            for word in re.findall(include_pattern, text_block):
                if word not in vocabulary and \
                        all([re.search(p, word) is None for p in RULES['exclude']['patterns']]):
                    i_line = len(text[:match.start() + 1].splitlines()) - 1
                    unknown_words[word]['files'][str(os.path.relpath(path, 'Spellcheck'))[3:]].append(i_line + 1)

        try:
            file = open(path, 'r')
            text = file.read()
            file_ext = path.suffix

            if file_ext in COMMENT_PATTERNS:
                file_pattern = COMMENT_PATTERNS[file_ext]
                for match in re.finditer(file_pattern, text):
                    block = text[match.start():match.end()]
                    _search_block(block)
            else:
                _search_block(text)

            file.close()

        except UnicodeDecodeError:
            errors.append(f'Could not read file: {path}')

    with IncrementalBar('Searching files', max=len(file_list)) as bar:
        for path in file_list:
            _search_file()
            bar.next()
        bar.finish()

    for error in errors:
        print(error)

    # TextBlob unknown words
    if len(unknown_words) > 0:
        os.system('python -m textblob.download_corpora')
        with IncrementalBar('Checking words', max=len(unknown_words)) as bar:
            for word in list(unknown_words.keys()):
                # TextBlob needs regular or lowercase word
                if word.isupper():
                    check_word = word.casefold()
                else:
                    check_word = word
                blob_words = Word(check_word).spellcheck()
                bwords_sorted = [[check_word, 0]]
                for w, c in blob_words:
                    c = round(c * 100, 2)
                    if w == word:
                        bwords_sorted[0][1] = c
                    else:
                        bwords_sorted.append([w, c])
                if RULES['auto_add_matches'] and bwords_sorted[0][1] == 100:
                    vocabulary.add(word)
                    unknown_words.pop(word)
                else:
                    unknown_words[word]['blob'] = bwords_sorted
                bar.next()
            bar.finish()
    return unknown_words


from pathlib import Path
import os
from collections import defaultdict
import json
import re
from typing import Union
from progress.bar import Bar
from textblob import Word
import pickle
from bz2 import BZ2File


PATH_VOCAB = 'tools/spellcheck/vocabulary.json'  # TODO: Path is required for PyGitHub
SP_DIR = os.path.dirname(os.path.realpath(__file__)) + '/'
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


RULES = load_json(Path(SP_DIR + 'rules.json'))
COMMENT_TYPES = load_json(Path(SP_DIR + 'comment_types.json'))
vocabulary = set(load_json(Path(SP_DIR + 'vocabulary.json')))


def write_json(obj, path: Path):
    """
    Write dict to .json file

    :param obj: Dictionary to save
    :param path: str of file path
    """
    with open(path, 'w') as file:
        json.dump(obj, file, indent=4)


def build_file_list(files_filter: Union[set, bool] = False) -> set:
    """
    Build list of all files to search, based on the rules provided.

    :param files_filter: Filter for included_files
    :return: set of included files
    """
    incl_ext = set(RULES['include']['extensions'])
    incl = set()
    excl_ext = set(RULES['exclude']['extensions'])
    for included in RULES['include']['locations']:
        for path in Path(included).rglob("*"):
            if path.is_file():
                file_ext = path.suffix
                if (len(incl_ext) > 0 and file_ext in incl_ext) or len(incl_ext) == 0:
                    if file_ext not in excl_ext \
                            and not any([path.is_relative_to(excl_loc) for excl_loc in RULES['exclude']['locations']]):
                        incl.add(path)
    file_list = files_filter.intersection(incl) if files_filter else incl
    return file_list


def get_words(files_filter: Union[set, bool] = False) -> defaultdict:
    """
    Find all valid words from all files defined by rules.json.

    :param files_filter: Filter for included_files
    :return: All valid words from all included files
    """
    unknown_words = defaultdict(lambda: {'error': '', 'action': {'replacement': '', 'vocabulary': ''},
                                         'blob': [], 'files': defaultdict(list)})
    errors = []
    file_list = build_file_list(files_filter)
    pattern = RULES['include']['pattern']

    def _search_file():

        def _search_block():
            for word in re.findall(pattern, txt_block):
                if word not in vocabulary and \
                        all([re.search(p, word) is None for p in RULES['exclude']['patterns']]):
                    unknown_words[word]['files'][str(os.path.relpath(path, 'Spellcheck'))[3:]].append(i_line + 1)

        try:
            file = open(path, 'r')
            txt = file.read()
            file_ext = path.suffix
            block_comment = ''
            for i_line, line in enumerate(txt.split('\n')):
                line_comment = ''
                for i_block, txt_block in enumerate(line.split()):
                    if file_ext in COMMENT_TYPES:
                        # Block comment end
                        if block_comment != '' and txt_block.endswith(block_comment):
                            txt_block = txt_block.strip(block_comment)
                            block_comment = ''
                            _search_block()
                        # Line comment option_start
                        if block_comment == '':
                            for lc in COMMENT_TYPES[file_ext]['line']:
                                i_lc = txt_block.find(lc)
                                if i_lc != -1:
                                    line_comment = lc
                                    txt_block = txt_block.strip(lc)
                                    break
                        # Block comment option_start
                        if i_block == 0 and line_comment == '':
                            for ext in COMMENT_TYPES[file_ext]['block']:
                                bc_start, bc_end = ext.split()
                                i_bc = txt_block.find(bc_start)
                                if 1 >= i_bc > -1:
                                    block_comment = bc_end
                                    txt_block = txt_block[i_bc:].strip(bc_start)
                                    break
                    if block_comment != '' or line_comment != '' or file_ext not in COMMENT_TYPES:
                        _search_block()

            file.close()
        except UnicodeDecodeError:  # TODO: Add error handling
            errors.append(f'{UnicodeDecodeError}: {path}')

    with Bar('Searching files..', max=len(file_list)) as bar:
        for path in file_list:
            _search_file()
            bar.next()
    if len(unknown_words) > 0:
        os.system('python -m textblob.download_corpora')
        for word in list(unknown_words.keys()):
            blob_words = Word(word).spellcheck()
            bwords_sorted = [[word, 0]]
            for w, c in blob_words:
                c = round(c * 100, 2)
                if w == word:
                    bwords_sorted[0][1] = c
                else:
                    bwords_sorted.append([w, c])
            unknown_words[word]['blob'] = bwords_sorted
    return unknown_words


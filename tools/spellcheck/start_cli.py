from spellcheck import *
import time
from PyInquirer import prompt, Separator
import fileinput


navigation = {
    'Main Menu': {
        'type': 'list',
        'name': 'main_menu',
        'message': 'Main Menu',
        'default': 'Settings',
        'choices': [
            'Search',
            {
                'name': 'Load & Update',
                'disabled': 'unknown_words.json missing',
                'value': 'Load & Update'
            },
            {
                'name': 'Start',
                'disabled': 'Use "Search" to find unknown words',
                'value': 'Start'
            },
            Separator(),
            'Save',
            {
                'name': 'Process All',
                'disabled': 'Nothing to process',
                'value': 'Process All'
            },
            'Exit'
        ]
    },
    'Start': {
        'type': 'list',
        'name': 'start',
        'message': '',
        'choices': [
            'Add to Vocabulary',
            'Replace',
            {
                'name': 'Skip',
                'disabled': '',
                'value': 'Skip'
            },
            'Jump',
            Separator(),
            'Main Menu',
            'Exit'
        ]
    },
    'Replace': {
        'type': 'input',
        'name': 'replace_input',
        'message': ''
    },
    'Replace all': {
        'type': 'confirm',
        'name': 'replace_all',
        'message': ''
    },
    'Vocabulary': {
        'type': 'input',
        'name': 'vocabulary',
        'message': ''
    },
    'Jump': {
        'type': 'list',
        'name': 'jump',
        'message': 'Jump to..',
        'choices': []
    },
    'Process All': {
        'type': 'confirm',
        'name': 'process',
        'message': 'Process all saved actions?'
    },
    'Exit': {
        'type': 'confirm',
        'name': 'exit',
        'message': 'Are you sure you want to exit the CLI? Progress will only be saved by using "Save" before Exit.'
    }
}

global vocabulary, reference, flat_vocab, vocab_keys, unknown_words
curr_id = 0
curr_word = ''


def add_to_vocab() -> str:
    """
    Interface for adding curr_word to vocabulary

    :return: Next menu branch: "Start"
    """
    print(vocab_keys)
    replacement = unknown_words[curr_word]['action']['replacement']
    navigation['Vocabulary']['message'] = f'Put "{curr_word if replacement == "" else replacement}" in..'
    category = list(prompt(navigation['Vocabulary']).values())[0]
    try:
        if 1 <= int(category) < len(reference):
            unknown_words[curr_word]['action']['vocab_index'] = category
            clear()
        else:
            raise ValueError
    except ValueError:
        print(f'Incorrect value, must be 1-{len(reference)}.')
        time.sleep(2)
        clear()
        add_to_vocab()
    return 'Start'


def replace():
    """
    Interface for replacing the current word

    :return: Next menu branch: "Start"
    """
    global curr_id

    for path, n_lines in unknown_words[curr_word]['files'].items():
        with open(path, 'r') as file:
            lines = file.readlines()
            for n in n_lines:
                for i in range(max(0, n-3), min(n+2, len(lines))):
                    print('[{x}]\t{y}'.format(x=(('\033[1m' + str(i) + '\033[0m') if i == n else i), y=lines[i]))
    navigation['Replace']['message'] = f'Replace "{curr_word}" with..'
    replacement = list(prompt(navigation['Replace']).values())[0]
    clear()
    navigation['Replace all']['message'] = f'Replace "{curr_word}" with "{replacement}"?'
    confirm = list(prompt(navigation['Replace all']).values())[0]
    clear()
    if confirm:
        unknown_words[curr_word]['action']['replacement'] = replacement
        if replacement not in flat_vocab:
            add_to_vocab()
        else:
            for i, words in reference.items():
                if replacement in set(words):
                    unknown_words[curr_word]['action']['vocab_index'] = i+1
    else:
        curr_id -= 1
    return 'Start'


def update_word_list():
    word_list = []
    for i, word in enumerate(list(unknown_words.keys())):
        replacement, vocab_index = unknown_words[word]['action'].values()
        jump_line = f'[{i + 1}] {word}'
        if replacement != '':
            navigation['Main Menu']['choices'][5]['disabled'] = ''
            jump_line += f' -> replace with: "{replacement}"'
        if vocab_index != '':
            navigation['Main Menu']['choices'][5]['disabled'] = ''
            jump_line += f' -> put in vocabulary index: {vocab_index}'
        word_list.append(jump_line)
    navigation['Jump']['choices'] = word_list


def process_actions():
    """
    Process all stored actions in unknown_words
    """
    for word in list(unknown_words.keys()):
        properties = unknown_words[word]
        replacement, vocab_index = properties['action'].values()
        if replacement != '':
            for path, lines in properties['files'].items():
                for n_line, line in enumerate(fileinput.input(path, inplace=True)):
                    if n_line+1 in lines:
                        print(line.replace(word, replacement), end='')
                    else:
                        print(line, end='')
        if vocab_index != '':
            reference[int(vocab_index)-1].append(word if replacement == '' else replacement)
            reference[int(vocab_index)-1].sort(key=str.casefold)
            unknown_words.pop(word)
    if len(unknown_words) == 0:
        os.remove(PATH_UNKNOWN_WORDS)
    print(unknown_words)
    write_json(VOCABULARY, Path(SP_DIR + 'vocabulary.json'))
    write_json(unknown_words, PATH_UNKNOWN_WORDS)
    print('All words processed, successfully!')


def update_word_header():
    header = f'{curr_id + 1}/{len(unknown_words)}: {curr_word}'
    replacement, vocab_index = unknown_words[curr_word]['action'].values()
    if replacement != '' or vocab_index != '':
        header += ' *'
    navigation['Start']['message'] = header


def start_cli():
    """
    Interface initiation
    """
    global curr_id, curr_word, unknown_words
    clear()

    branch = 'Main Menu'
    while True:
        clear()

        branch = list(prompt(navigation[branch]).values())[0]

        # "Main Menu" and Menu choices
        if branch == 'Main Menu' or branch == 'Start':
            pass
        elif branch == 'Search':
            unknown_words = get_words()
            print(f'Found {len(unknown_words)} unknown words!')
            time.sleep(2)
            navigation['Main Menu']['choices'][2]['disabled'] = ''
            branch = 'Main Menu'
        elif branch == 'Load & Update':

            unknown_words = get_words()
            old_unknown_words = load_json(PATH_UNKNOWN_WORDS)
            for word, properties in old_unknown_words.items():
                if word in unknown_words:
                    unknown_words[word]['action'] = properties['action']

            print(f'Loaded old unknown words file!')
            time.sleep(2)
            navigation['Main Menu']['choices'][2]['disabled'] = ''
            branch = 'Main Menu'
        elif branch == 'Exit':
            confirm = list(prompt(navigation[branch]).values())[0]
            if confirm:
                exit()
            else:
                branch = 'Main Menu'
        elif branch == 'Save':
            write_json(unknown_words, Path('unknown_words.json'))
            print('Progress saved!')
            branch = 'Main Menu'
            time.sleep(2)
        elif branch == 'Process All':
            clear()
            confirm = list(prompt(navigation['Process All']).values())[0]
            if confirm:
                process_actions()
            branch = 'Main Menu'
            clear()

        else:  # special choices inside "Start"
            if branch == 'Add to Vocabulary':
                branch = add_to_vocab()
            elif branch == 'Replace':
                branch = replace()
            elif branch == 'Skip':
                branch = 'Start'
            elif branch == 'Jump':
                curr_id = int(list(prompt(navigation['Jump']).values())[0].split()[0][1:-1]) - 2
                clear()
                branch = 'Start'

            curr_id += 1
            curr_id = min([curr_id, len(unknown_words)-1])
        navigation['Start']['choices'][2]['disabled'] = '' if curr_id+1 != len(unknown_words) else 'end'
        curr_word = list(unknown_words.keys())[curr_id]
        update_word_header()
        update_word_list()


def main():
    global vocabulary, reference, flat_vocab, vocab_keys

    # Load all
    reference = set_ref()
    flat_vocab = flatten_vocab()
    vocab_keys = get_vocab_keys('::', '......')

    if os.path.exists(PATH_UNKNOWN_WORDS):
        navigation['Main Menu']['choices'][1]['disabled'] = ''

    # bring up user prompt
    start_cli()


if __name__ == '__main__':
    main()

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
        'name': 'option_start',
        'message': '',
        'choices': []
    },
    'Continue': {
        'type': 'input',
        'name': 'show_occurrences',
        'message': 'Press any key to continue..'
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

global unknown_words
curr_id = 0
curr_word = ''


def show_occurrences():
    for path, n_lines in unknown_words[curr_word]['files'].items():
        print('\033[1m' + str(path) + '\033[0m')
        with open(path, 'r') as file:
            lines = file.readlines()
            for n in n_lines:
                for i in range(max(0, n-5), min(n+4, len(lines))):
                    print('[{x}]\t{y}'.format(x=(('\033[1m' + str(i) + '\033[0m') if i == n-1 else i), y=lines[i]))

    prompt(navigation['Continue'])


def add_action(replacement):
    unknown_words[curr_word]['action']['replacement'] = replacement
    if replacement not in vocabulary:
        unknown_words[curr_word]['action']['vocabulary'] = replacement


def replace():
    global curr_id
    navigation['Replace']['message'] = f'Replace "{curr_word}" with..'
    replacement = list(prompt(navigation['Replace']).values())[0]
    clear()
    navigation['Replace all']['message'] = f'Replace "{curr_word}" with "{replacement}"?'
    confirm = list(prompt(navigation['Replace all']).values())[0]
    clear()
    if confirm:
        add_action(replacement)
    else:
        curr_id -= 1


def update_word_list():
    word_list = []
    for i, word in enumerate(list(unknown_words.keys())):
        replacement, vocab_word = unknown_words[word]['action'].values()
        jump_line = f'[{i + 1}] {word}'
        if replacement != '':
            navigation['Main Menu']['choices'][5]['disabled'] = ''
            jump_line += f' -> replace with: "{replacement}"'
        if vocab_word != '':
            navigation['Main Menu']['choices'][5]['disabled'] = ''
            jump_line += f' -> put in vocabulary'
        word_list.append(jump_line)
    navigation['Jump']['choices'] = word_list


def update_word_header():
    header = f'{curr_id + 1}/{len(unknown_words)}: {curr_word}'
    replacement, vocab_index = unknown_words[curr_word]['action'].values()
    if replacement != '' or vocab_index != '':
        header += ' *'
    navigation['Start']['message'] = header


def update_replace_list():
    replace_list = [f'"{word}": {certainty}%' for word, certainty in unknown_words[curr_word]['blob']]
    replace_list.extend([Separator(), 'Custom', 'Add to Vocabulary', Separator(), 'Show occurrences',
                         {'name': 'Skip', 'disabled': '', 'value': 'Skip'}, 'Jump', 'Main Menu', 'Exit'])
    navigation['Start']['choices'] = replace_list


def process_actions():
    """
    Process all stored actions in unknown_words
    """
    global curr_id

    for i_word, word in enumerate(list(unknown_words.keys())):
        properties = unknown_words[word]
        replacement, vocab_word = properties['action'].values()
        if replacement != '':
            for path, lines in properties['files'].items():
                for n_line, line in enumerate(fileinput.input(path, inplace=True)):
                    if n_line+1 in lines:
                        print(line.replace(word, replacement), end='')
                    else:
                        print(line, end='')
        if vocab_word != '':
            vocabulary.add(vocab_word)
            unknown_words.pop(word)
            if curr_id >= i_word:
                curr_id -= 1
    if len(unknown_words) == 0:
        os.remove(PATH_UNKNOWN_WORDS)
    new_vocab = list(vocabulary)
    new_vocab.sort()
    write_json(new_vocab, Path(SP_DIR + 'vocabulary.json'))
    print('All words processed, successfully!')


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

        # Choices
        if branch == 'Start':
            pass
        elif branch == 'Main Menu':
            pass
        elif branch == 'Search':
            unknown_words = get_words()
            print(f'Found {len(unknown_words)} unknown words!')
            navigation['Main Menu']['choices'][2]['disabled'] = ''
            branch = 'Main Menu'
            prompt(navigation['Continue'])
        elif branch == 'Load & Update':

            unknown_words = get_words()
            old_unknown_words = load_json(PATH_UNKNOWN_WORDS)
            for word, properties in old_unknown_words.items():
                if word in unknown_words:
                    unknown_words[word]['action'] = properties['action']

            print(f'Loaded old unknown words file!')
            prompt(navigation['Continue'])
            navigation['Main Menu']['choices'][2]['disabled'] = ''
            branch = 'Main Menu'
        elif branch == 'Exit':
            confirm = list(prompt(navigation[branch]).values())[0]
            if confirm:
                exit()
            else:
                branch = 'Main Menu'
        elif branch == 'Save':
            write_json(unknown_words, Path(PATH_UNKNOWN_WORDS))
            print('Progress saved!')
            branch = 'Main Menu'
            prompt(navigation['Continue'])
        elif branch == 'Process All':
            clear()
            confirm = list(prompt(navigation['Process All']).values())[0]
            if confirm:
                process_actions()
                prompt(navigation['Continue'])
            branch = 'Main Menu'
            clear()

        else:  # special choices inside "Start"
            if branch == 'Skip':
                pass
            elif branch == 'Jump':
                curr_id = int(list(prompt(navigation['Jump']).values())[0].split()[0][1:-1]) - 2
                clear()
            elif branch == 'Show occurrences':
                show_occurrences()
                curr_id -= 1
            elif branch == 'Custom':
                replace()
            elif branch == 'Add to Vocabulary':
                unknown_words[curr_word]['action']['vocabulary'] = curr_word
            else:
                replacement = branch[1:].split('"')[0]
                add_action(replacement)
            branch = 'Start'
            curr_id += 1
            curr_id = min([curr_id, len(unknown_words)-1])
        curr_word = list(unknown_words.keys())[curr_id]
        update_word_header()
        update_word_list()
        update_replace_list()
        navigation['Start']['choices'][-4]['disabled'] = '' if curr_id + 1 != len(unknown_words) else 'end'


def main():
    if os.path.exists(PATH_UNKNOWN_WORDS):
        navigation['Main Menu']['choices'][1]['disabled'] = ''

    # bring up user prompt
    start_cli()


if __name__ == '__main__':
    main()

import os.path
from utils import *
from PyInquirer import prompt, Separator
import fileinput

# How many surrounding lines should be printed
SCAN_WINDOW = 10

navigation = {
    'Main Menu': {
        'type': 'list',
        'name': 'main_menu',
        'message': 'Main Menu',
        'default': 'Settings',
        'choices': [
            'Search',
            {
                'name': 'Load',
                'disabled': 'unknown_words.json missing',
                'value': 'Load'
            },
            {
                'name': 'Update',
                'disabled': 'Nothing to update',
                'value': 'Update'
            },
            {
                'name': 'Start',
                'disabled': 'No words to start',
                'value': 'Start'
            },
            Separator(),
            {
                'name': 'Save',
                'disabled': 'Nothing to save',
                'value': 'Save'
            },
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
    """
    Print the surrounding file content of the word
    """
    for path, n_lines in unknown_words[curr_word]['files'].items():
        print('\033[1m' + str(path) + '\033[0m')
        with open(path, 'r') as file:
            lines = file.readlines()
            for n in n_lines:
                for i in range(max(0, n-SCAN_WINDOW), min(n+SCAN_WINDOW, len(lines))):
                    if i == n - 1:
                        prefix = '->\t'
                        number = ('\033[1m' + str(i) + '\033[0m')
                    else:
                        prefix = '\t'
                        number = str(i)
                    print('{a} [{b}]\t{c}'.format(a=prefix, b=number, c=lines[i]))

    prompt(navigation['Continue'])


def assign_action(replacement):
    """
    Assign replacement to word
    """
    unknown_words[curr_word]['action']['replacement'] = replacement
    if replacement not in vocabulary:
        unknown_words[curr_word]['action']['vocabulary'] = replacement


def replace():
    """
    Custom word replacement
    """
    global curr_id
    navigation['Replace']['message'] = f'Replace "{curr_word}" with..'
    replacement = list(prompt(navigation['Replace']).values())[0]
    clear()
    navigation['Replace all']['message'] = f'Replace "{curr_word}" with "{replacement}"?'
    confirm = list(prompt(navigation['Replace all']).values())[0]
    clear()
    if confirm:
        assign_action(replacement)
    else:
        curr_id -= 1


def update_word_list():
    """
    Update the word list when using "Jump"
    """
    word_list = []
    for i, word in enumerate(list(unknown_words.keys())):
        replacement, vocab_word = unknown_words[word]['action'].values()
        jump_line = f'[{i + 1}] {word}'
        if replacement != '':
            navigation['Main Menu']['choices'][6]['disabled'] = ''
            jump_line += f' -> replace with: "{replacement}"'
        if vocab_word != '':
            navigation['Main Menu']['choices'][6]['disabled'] = ''
            jump_line += f' -> put in vocabulary'
        word_list.append(jump_line)
    navigation['Jump']['choices'] = word_list


def update_word_header():
    """
    Update the word header shown in "Start"
    """
    header = f'{curr_id + 1}/{len(unknown_words)}: {curr_word}'
    replacement, vocab_index = unknown_words[curr_word]['action'].values()
    if replacement != '' or vocab_index != '':
        header += ' *'
    navigation['Start']['message'] = header


def update_start():
    """
    Update "Start" choices
    """
    bwords = unknown_words[curr_word]["blob"]
    replace_list = [f'Add to Vocabulary ({bwords[0][1]}%)']
    if len(bwords) > 1:
        replace_list.extend([f'Replace with "{word}" ({certainty}%)' for word, certainty in bwords[1:]])
    replace_list.extend(['Replace with custom', Separator(), 'Show occurrences',
                         {'name': 'Skip', 'disabled': '', 'value': 'Skip'}, 'Jump', 'Main Menu', 'Exit'])
    navigation['Start']['choices'] = replace_list


def disable_choices():
    """
    Disable all Main Menu choices depending on unknown words
    """
    navigation['Main Menu']['choices'][1]['disabled'] = 'unknown_words.json missing'
    navigation['Main Menu']['choices'][2]['disabled'] = 'Nothing to update'
    navigation['Main Menu']['choices'][3]['disabled'] = 'No words to start'
    navigation['Main Menu']['choices'][5]['disabled'] = 'Nothing to save'
    navigation['Main Menu']['choices'][5]['disabled'] = 'Nothing to process'


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
            if curr_id > i_word:
                curr_id -= 1
    if len(unknown_words) == 0 and os.path.exists(PATH_UNKNOWN_WORDS):
        os.remove(PATH_UNKNOWN_WORDS)
        disable_choices()
    save_vocab()
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
            unknown_words = get_words('../../')
            save_vocab()
            print(f'Found {len(unknown_words)} unknown words!')
            branch = 'Main Menu'
            prompt(navigation['Continue'])
            if len(unknown_words) == 0:
                continue
            navigation['Main Menu']['choices'][2]['disabled'] = ''
            navigation['Main Menu']['choices'][3]['disabled'] = ''
            navigation['Main Menu']['choices'][5]['disabled'] = ''

        elif branch == 'Load':
            unknown_words = load_json(PATH_UNKNOWN_WORDS)
            print(f'Loaded {len(unknown_words)} unknown words!')
            prompt(navigation['Continue'])
            branch = 'Main Menu'
            navigation['Main Menu']['choices'][2]['disabled'] = ''
            navigation['Main Menu']['choices'][3]['disabled'] = ''
            navigation['Main Menu']['choices'][5]['disabled'] = ''

        elif branch == 'Update':
            new_unknown_words = get_words('../../')
            save_vocab()
            for word, properties in unknown_words.items():
                if word in new_unknown_words:
                    unknown_words[word]['action'] = properties['action']

            print(f'Updated unknown words!')
            prompt(navigation['Continue'])
            navigation['Main Menu']['choices'][3]['disabled'] = ''
            branch = 'Main Menu'
        elif branch == 'Exit':
            confirm = list(prompt(navigation[branch]).values())[0]
            if confirm:
                exit()
            else:
                branch = 'Main Menu'
        elif branch == 'Save':
            save_json(unknown_words, Path(PATH_UNKNOWN_WORDS))
            print('Progress saved!')
            branch = 'Main Menu'
            navigation['Main Menu']['choices'][1]['disabled'] = ''
            prompt(navigation['Continue'])
        elif branch == 'Process All':
            clear()
            confirm = list(prompt(navigation['Process All']).values())[0]
            if confirm:
                process_actions()
                prompt(navigation['Continue'])
            branch = 'Main Menu'
            clear()
            if len(unknown_words) == 0:
                print('No unknown words left.')
                disable_choices()
                prompt(navigation['Continue'])
                continue

        else:  # "Start" choices
            if branch == 'Skip':
                pass
            elif branch == 'Jump':
                curr_id = int(list(prompt(navigation['Jump']).values())[0].split()[0][1:-1]) - 2
                clear()
            elif branch == 'Show occurrences':
                show_occurrences()
                curr_id -= 1
            elif branch.startswith('Replace with custom'):
                replace()
            elif branch.startswith('Add to Vocabulary'):
                unknown_words[curr_word]['action']['vocabulary'] = curr_word
            else:
                replacement = branch[1:].split('"')[1]
                assign_action(replacement)
            curr_id += 1
            if curr_id < len(unknown_words):
                branch = 'Start'
            else:
                branch = 'Main Menu'
            curr_id = min([curr_id, len(unknown_words)-1])
        curr_word = list(unknown_words.keys())[curr_id]
        update_word_header()
        update_word_list()
        update_start()
        navigation['Start']['choices'][-4]['disabled'] = '' if curr_id + 1 != len(unknown_words) else 'end'


def main():
    if os.path.exists(PATH_UNKNOWN_WORDS):
        navigation['Main Menu']['choices'][1]['disabled'] = ''

    # bring up user prompt
    start_cli()


if __name__ == '__main__':
    main()

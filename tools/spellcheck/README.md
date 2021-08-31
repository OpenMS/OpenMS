# Spellcheck


## Set-Up

### Rules
To work properly, `rules.json` must be set up correctly.
- `include`
  - `locations` Directories to apply search to (Relative to the working directory)
  - `pattern` Pattern to classify phrase as word
  - `extensions` (Optional) file types to include in search, e.g. `".py"`
- `exclude` 
  - `locations` (Optional) exclude files or directories (Relative to the working directory)
  - `patterns` (Optional) exclude patterns
  - `extensions` (Optional) exclude file types

### Comment Types
For each file type, where only comments should be searched for words, `line` and `block` comment types must be provided
in `comment_types.json`.  
Example: 
```
".py": {
        "line":  
        [
            "#"
        ],
        "block": [
            "''' '''",
            "\"\"\" \"\"\""
        ]
    }
```


## Command Line Interface (CLI)

### Steps
1. Install Python >= 3.8.x
2. `pip install -r requirements.txt`
3. Run `python start_cli.py` from command window
4. Use the CLI to process all found unknown words

### Additional Info
No decision is executed directly, instead all actions are stored until `Process All` is used.
- `Search` triggers a search for unknown words over all files specified by the rules in `rules.json`.
- `Load & Update` loads an already saved file of unknown words, but will also trigger a new search. 
  Previously stored actions are then assigned to identical words.
  This is because its occurrences may change, or the word may have been deleted. 
- `Start` will start with the first detected word and continue with the remaining unknown words.
  Already assigned words are marked with a `*` at the end of the header. The assignment can be viewed by using `Jump` to
  show all found words and their assigned actions.
- `Process All` will trigger all stored decisions to be executed. 
  This option is activated after at least one action is assigned to at least one unknown word.
- `Save` will write the current progress to `unknown_words.json` in the working directory.
- `Exit` will close the program, progress will not be saved automatically!


## GitHub
Unknown words are searched automatically by using GitHub actions workflows, every time a new commit is pushed. 
If new words are found, the GitHub bot will open an issue or update an existing one on the respective fork or repository.
By default the action only searches updated files, a full search can be triggered by clicking on `Actions` and 
`spellcheck - show unknown words` in the GitHub repository, where the search should be conducted.  
The open issue allows to directly assign actions to identified words, by editing the issue and replacing the whitespace 
between \` \` in the ` `-code field.  
**Editing and saving the issue will trigger all actions to be processed inside the repository!**  
The workflow will always search and process in the edited branch.


## Author
**Matteo Pilz**

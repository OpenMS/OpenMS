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
- `auto_add_matches` <`true`/`false`> Checked words with a 100% certainty to be correct are automatically placed in the
vocabulary

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
3. Run `python start_cli.py` from command window inside `<Project>/.github/spellcheck`
4. Use the CLI to process all found unknown words

### Additional Info
No decision is executed directly, instead all actions are stored until `Process All` is used.
- `Search` triggers a search for unknown words over all files specified by the rules in `rules.json`.
- `Load` loads an already saved file of unknown words.  
   **IMPORTANT:** If files have been changed after saving, it is necessary to use `Update` after `Load`
- `Update` will trigger a new search and assign previously stored actions to the updated words.  
  Updating is necessary, as word occurrences may change, or it may have been deleted. 
- `Start` will start with the first detected word and continue with the remaining unknown words.
  Already assigned words are marked with a `*` at the end of the header. The assignment can be viewed by using `Jump`.
- `Process All` will trigger all stored decisions to be executed. 
  This option is enabled after at least one action is assigned to at least one unknown word.
- `Save` will write the current progress to `unknown_words.json` in the working directory.
- `Exit` will close the program, progress will not be saved automatically!


## GitHub

### General
Unknown words are searched automatically by using GitHub actions workflows, every time a new commit is pushed. 
If new words are found, the GitHub bot will open an issue or update an existing one on the respective fork or repository.
By default, the action only searches edited files. A full search can be initiated manually by clicking on `Actions` and 
`Spellcheck - Full Search` in the GitHub repository. Additionally, the correct branch needs to be selected.

### Assigning Actions
The open issue allows to directly assign actions to identified words. Using the checkboxes, words can be either replaced
or added to the vocabulary *(Check only one option per word)*. To replace with a custom word, the issue needs to be edited
and the underscore in the line: `Replace "<word>" with "_" (custom)` replaced with the correct word.

### Processing
To trigger the processing of all decisions, create a comment with the word "process" in the issue. Processing in GitHub
will edit the respective files and the vocabulary directly in the repository.


## Author
**Matteo Pilz**

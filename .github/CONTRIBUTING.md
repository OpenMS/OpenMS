OpenMS Documentation Contributing Guidelines
============================================

Hi there!

Thank you for thinking about enhancing OpenMS Documentation further.

Please feel free to:

1. Create a bug report or feature request in OpenMS Documentation, [here](https://github.com/OpenMS/OpenMS-docs/issues).
2. Create a pull request in [OpenMS Documentation](https://github.com/OpenMS/OpenMS-docs) with the change you're proposing.

For any questions, drop us a ping in [OpenMS Gitter](https://gitter.im/OpenMS/OpenMS).

## Create a Pull Request(PR)

1. Fork this repository.
2. Add the change in your fork.
3. Create a pull request to [OpenMS/OpenMS-Docs](https://github.com/OpenMS/OpenMS-docs) `develop` branch.
4. Make sure you attach a few screenshots as to how your change looks like.
5. If this change belongs to [OpenMS API reference](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/index.html),
   please create a pull request [here](https://github.com/OpenMS/OpenMS/tree/develop/doc).

## Documentation content writing conventions

1. Please prefer active voice.
2. Restrict the line length to 120 characters, including space(s).
3. Title and H1 heading should be in `Title Case`. Follow `Sentence case`, until otherwise stated.
4. Write as OpenMS is writing for itself as an object and subject.
5. Use American English.
6. Use `backtick`(s) for formatting code, library name, files, and path.
7. Use **bold** for product name, object name, an independent entity.
8. Use **bold** for menu title in an application.
9. Link to glossary terms using {term}\`this is a glossary term\`.
10. OpenMS documentation uses following Admonitions
    - Hint
	- Important
	- Note
	- Warning
	- Tip
	- See Also

	Example of these are present in the documentation, please follow them.
11. OpenMS documentation features admonitions for links to videos. If you wish to add a link to a video, use following html code:

	<div class="admonition video">
		<p class="admonition-title">**Video**</p>
		<!-- Insert text -->
	</div>
12. Always specify lexers for code blocks.
13. Format keyboard strokes using `<kbd>qwerty-keyboard-button</kbd>`, as an example, please see [this](../docs/tutorials/TOPP/hotkeys-table.md) file.
14. Be nice, polite, and respectful.


### Naming files

1. The title of the page should be the `name-of-the-file.md`.
2. Prefer writing in markdown, with an `.md` file extension; both reStructureText and markdown is supported in OpenMS Documentation.

### Images and figures

For images and figures:

1. Add a screenshot of the window.
2. In tutorial, align images in center. Other instructions, should have aligment to left.
3. Please set the size to `500px` of images added in step-by-step guides or instructions.

## OpenMS documentation contributors

Thank you for your contribution!

Finally, please add your name below:

1. OpenMS Team
2. Michael R. Crusoe

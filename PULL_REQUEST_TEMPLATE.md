## Description

<!-- Please include a summary of the change and which issue is fixed here. -->

## Checklist
- [ ] Make sure that you are listed in the AUTHORS file
- [ ] Add relevant changes and new features to the CHANGELOG file
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] New and existing unit tests pass locally with my changes
- [ ] Updated or added python bindings for changed or new classes (Tick if no updates were necessary.)

### How can I get additional information on failed tests during CI
<details>
  <summary>Click to expand</summary>
If your PR is failing you can check out

- The details of the action statuses at the end of the PR or the "Checks" tab.
- http://cdash.seqan.de/index.php?project=OpenMS and look for your PR. Use the "Show filters" capability on the top right to search for your PR number.
  If you click in the column that lists the failed tests you will get detailed error messages.

</details>

### Advanced commands (admins / reviewer only)
<details>
  <summary>Click to expand</summary>
  
- `/reformat` (experimental) applies the clang-format style changes as additional commit. Note: your branch must have a different name (e.g., yourrepo:feature/XYZ) than the receiving branch (e.g., OpenMS:develop). Otherwise, reformat fails to push.
- setting the label "NoJenkins" will skip tests for this PR on jenkins (saves resources e.g., on edits that do not affect tests)
- commenting with `rebuild jenkins` will retrigger Jenkins-based CI builds
  
</details>

---
:warning: Note: Once you opened a PR try to minimize the number of *pushes* to it as every push will trigger CI (automated builds and test) and is rather heavy on our infrastructure (e.g., if several pushes per day are performed).

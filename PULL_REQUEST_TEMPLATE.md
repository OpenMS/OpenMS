# Description

Please include a summary of the change and which issue is fixed.

# Checklist:
- [ ] Make sure that you are listed in the AUTHORS file
- [ ] Add relevant changes and new features to the CHANGELOG file
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] New and existing unit tests pass locally with my changes
- [ ] Updated or added python bindings for changed or new classes. (Tick if no updates were necessary.)

# How can I get additional information on failed tests during CI:
If your PR is failing you can check out 
- http://cdash.openms.de/index.php?project=OpenMS and look for your PR. If you click in the column that lists the failed tests you will get detailed error messages.
- Or click on the action: e.g., for clang-format linting

# Note:
- Once you opened a PR try to minimize the number of *pushes* to it as every push will trigger CI (automated builds and test) and is rather heavy on our infrastructure (e.g., if several pushes per day are performed).

# Advanced commands (admins / reviewer only):
- /rebase will try to rebase the PR on the current develop branch.
- /reformat (experimental) applies the clang-format style changes as additional commit
- setting the label NoJenkins will skip tests for this PR on jenkins (saves resources e.g., on edits that do not affect tests)

Pull Request Checklist
======================

Before opening a pull request, check the following:

1. **Does the code build?**
  Execute `make` (or your build system's equivalent, e.g., `cmake --build . --target ALL_BUILD --config Release` on 
  Windows).
2. **Do all tests pass?**
   To check if all tests have passed, execute `ctest`.
   If a test that is unrelated to your changes fails, check the [nightly builds](https://cdash.seqan.de/index.php?project=OpenMS)
   to see if the error is also in `develop`. If the error is in develop, [create a github issue](/manual/contribute.md#Write and Label GitHub Issues).
3. **Is the code documented?**
   Document all new classes, including their methods and parameters.
   It is also recommended to document non-public members and methods.
4. **Does the code introduce changes to the API?**
   If the code introduces changes to the API, make sure that the documentation is up-to-date and that the Python bindings
   (pyOpenMS) still work. For each change in the C++ API, make a change in the Python API wrapper via 
   the `pyOpenMS/pxds/` files.
5. **Have you completed regression testing?**
   Make sure that you include a test in the test suite for:
   - Public methods of a class
   - TOPP tools
   - Bug fixes

Make sure to:

- **Rebase before you open a pull request.**
  To include all recent changes, rebase your branch on `develop` before opening a pull request.
  If you pushed your branch to `origin` before rebasing, git will most likely tell you after the rebase that your
  local branch and the remote branch have diverged. If you are sure that the remote branch does not contain any local
  commits in the rebased version, you can safely push using `git push -f origin <branch-name>` to enforce overwrite. If
  not, contact your local git expert on how to get the changes into your local branch.

- **Capture similar changes in a single commit**
  Each commit should represent one logical unit. Consolidate multiple commits if they belong together or split single
  commits if they are unrelated. For example, committing code formatting together with a one-line fix makes it very hard
  to figure out what the fix was and which changes were inconsequential.

* **Create a pull request for a single feature or bug**
  If you have multiple features or fixes in a pull request, you might get asked to split your request and open multiple
  pull requests instead.

* **Describe what you have changed in your pull request.**
  When opening the pull request, give a detailed overview of what has changed and why. Include a clear rationale for the
  changes and add benchmark data if available. See [this request](https://github.com/bitly/dablooms/pull/19) for 
  an example.

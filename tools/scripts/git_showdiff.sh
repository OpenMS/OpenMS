
echo "Usage: git_showdiff.sh REV1 REV2"

git diff $1 $2 --name-only \
  | grep -v '^src/tests/topp/' \
  | grep -v '^share/OpenMS' \
  | grep -v '^src/tests/class_tests/openms/data/' \
  | grep -v '^cmake/modules/'  \
  | grep -v '^cmake/'  \
  | grep -v '^src/openswathalgo/thirdparty' \
  | grep -v '^src/openms/thirdparty/' \
  | grep -v '^tools'  \
  | xargs git diff --shortstat $1 $2 --



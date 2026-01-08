import os
import subprocess
import glob

# script that lists all test files in the topp folder that are:
# 1. tracked by git
# 2. not referenced in the topp tests (CMakeList.txt)
# after careful check these potentially can be removed

def list_files_in_directory(start_path):
    """ Returns a list of all files in a directory (no recurse). """
    files_list = []
    for file in os.listdir(start_path):
        full_path = os.path.join(start_path, file)
        if os.path.isfile(full_path):
            files_list.append(os.path.join(start_path, file))
    return files_list

def file_contains_text(file_path, text):
    """ Checks if the given text is in the file. """
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            return text in file.read()
    except UnicodeDecodeError:
        return False

def resolve_files(file_list):
    """ Resolves a list of file paths and wildcard patterns to actual file paths. """
    resolved_files = []
    for file_path in file_list:
        # Check if the file path contains a wildcard
        if '*' in file_path:
            # Expand the wildcard pattern to actual file paths
            expanded_files = glob.glob(file_path, recursive=True)
            resolved_files.extend(expanded_files)
        else:
            # If no wildcard, just add the file path as is
            if os.path.exists(file_path):
                resolved_files.append(file_path)
            else:
                print(f"Warning: The file '{file_path}' does not exist.")
    return resolved_files


def main(test_data_dir, source_files_to_check):
    test_files = list_files_in_directory(test_data_dir)
    source_files = resolve_files(source_files_to_check) # glob wildcards

    test_files_not_in_sources = []

    for test_file in test_files:
        test_file_name = os.path.basename(test_file)
        found = False
        for source_file in source_files:
            if file_contains_text(source_file, test_file_name):
                found = True
                break
        if not found: 
            test_files_not_in_sources.append(test_file)

    #print("Test files not found in any source file:")
    #for file in test_files_not_in_sources:
    #    print(file)

    print("Test files tracked by GIT but not found in any source file:")
    tracked = filter_tracked_files(test_files_not_in_sources) # print tracked files not referenced in sources
    for file in tracked:
        print(file)
    
def get_tracked_files():
    """Returns a set of tracked files in the current Git repository."""
    try:
        # Using git ls-files to list tracked files
        result = subprocess.run(['git', 'ls-files', test_data_directory], capture_output=True, text=True)
        result.check_returncode()  # Ensures that the git command didn't fail
        # Split the output by lines to get individual files
        tracked_files = set(result.stdout.splitlines())
        tracked_files = [os.path.basename(file) for file in tracked_files]
        return tracked_files
    except subprocess.CalledProcessError as e:
        print("Failed to run git command:", e)
        return set()

def filter_tracked_files(file_list):
    """Filters the given list of files to find which are tracked."""
    tracked_files = get_tracked_files()
    # Find the intersection of our file list with the tracked files
    return [file for file in file_list if os.path.basename(file) in tracked_files]

if __name__ == "__main__":
    # file names in test data will be checked for existance in source file
    test_data_directory = '../../src/tests/topp'
    source_files_to_check = [ '../../src/tests/topp/CMakeLists.txt',  '../../src/tests/class_tests/openms/source/*.cpp']
    main(test_data_directory, source_files_to_check)


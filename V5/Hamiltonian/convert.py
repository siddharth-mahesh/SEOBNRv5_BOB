import sys
import re

def read_file(file_path):
    """
    Reads a Python file and returns its content as a string.

    :param file_path: Path to the Python file to read.
    :return content: Content of the file as a string.
    :raises FileNotFoundError: If the specified file does not exist.
    """
    with open(file_path, 'r') as file:
        content = file.read()
    return content

def replace_mpf_with_sp_rational(content):
    """
    Replaces all occurrences of 'mpf(' with 'sp.Rational(' in the given content string.

    :param content: String content to modify.
    :return modified_content: Modified string with replacements made.
    """
    return content.replace('mpf(', 'sp.Rational(')

def replace_sp_rational_integers(content):
    """
    Replaces all occurrences of 'sp.Rational([integer].0[not an integer])' with '[integer]' in the given content string.

    :param content: String content to modify.
    :return modified_content: Modified string with replacements made.
    """
    pattern = re.compile(r'sp\.Rational\((\d+)\.0([^\d])\)')
    return pattern.sub(r'\1\2', content)

def replace_integers_with_fractional_parts(content):
    """
    Replaces all numbers matching the pattern '[integer].0[not an integer]' with '[integer]' and
    '[integer].5[not an integer]' with 'sp.Rational(2*[integer]+1, 2)' in the given content string.

    :param content: String content to modify.
    :return modified_content: Modified string with replacements made.
    """
    content = re.sub(r'(\d+)\.0([^\d])', r'\1\2', content)
    content = re.sub(r'(\d+)\.5([^\d])', r'sp.Rational(2*\1+1, 2)\2', content)
    return content

def write_to_file(content, original_file_path):
    """
    Writes the modified content to a new file with '-converted' appended to the original file name.

    :param content: String content to write to the file.
    :param original_file_path: Original file path to derive the new file name.
    :return: None
    """
    new_file_path = original_file_path.replace('.py', '_converted.py')
    with open(new_file_path, 'w') as file:
        file.write(content)

def main():
    """
    Main function to execute the script steps in order.

    :return: None
    :raises ValueError: If no file path is provided as command line argument.
    """
    if len(sys.argv) < 2:
        raise ValueError("No file path provided. Please specify a Python file to process.")

    file_path = sys.argv[1]

    # Step 1: Read the file
    content = read_file(file_path)

    # Step 2: Replace 'mpf(' with 'sp.Rational('
    content = replace_mpf_with_sp_rational(content)

    # Step 3: Replace 'sp.Rational([integer].0[not an integer])' with '[integer]'
    content = replace_sp_rational_integers(content)

    # Step 4: Replace '[integer].0[not an integer]' with '[integer]' and '[integer].5[not an integer]' with 'sp.Rational(2*[integer]+1, 2)'
    content = replace_integers_with_fractional_parts(content)

    # Step 5: Write the modified content to a new file
    write_to_file(content, file_path)

if __name__ == "__main__":
    main()

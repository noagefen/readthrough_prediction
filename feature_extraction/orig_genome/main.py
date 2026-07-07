import RNA
import os
import re
import tempfile


def process_RNA_sequences_v3(folder_path):
    """
    Process RNA sequences in all files in the specified folder.
    Each file is expected to be in FASTA format.
    The results are saved in a corresponding output file with '_res.fas' suffix.
    """

    def extract_min_value_and_index(output):
        """
        Extract the minimum value and its index from the RNALFold output.
        """
        min_value = float('inf')
        min_index = -1

        # Updated regex pattern
        pattern = re.compile(r'\(\s*(-?\d+\.\d+)\)\s+(\d+)$')

        for line in output.split('\n'):
            match = pattern.search(line)
            if match:
                value = float(match.group(1))
                index = int(match.group(2))  # Adjusted to the correct group number
                if value < min_value:
                    min_value = value
                    min_index = index

        return min_value, min_index

    for filename in os.listdir(folder_path):
        if filename.endswith(".fas"):
            with open(os.path.join(folder_path, filename), 'r') as file:
                lines = file.readlines()

            output_filename = filename.replace('.fas', '_res.fas')
            with open(os.path.join(folder_path,'output', output_filename), 'w') as output_file:
                sequence_name = ""
                sequence = ""

                for line in lines:
                    if line.startswith('>'):
                        # Process the previous sequence if it exists
                        if sequence:
                            # Creating a temporary file for RNALFold output
                            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
                                temp_file_name = temp_file.name

                                # Call the RNALFold function, passing the sequence and the temp file name
                                # Replace this with actual call to RNALFold(sequence, temp_file_name)
                                RNA.Lfold(sequence, 150, temp_file)

                                # RNALFold(sequence, temp_file_name)

                                # Read the output from the temporary file
                                temp_file.seek(0)
                                rnal_fold_output = temp_file.read()

                                min_value, min_index = extract_min_value_and_index(rnal_fold_output)
                                output_file.write(
                                    f">>{sequence_name}\nmax energy = {min_value} index = {min_index}\n\n")

                        sequence_name = line.strip()[1:]
                        sequence = ""
                    else:
                        sequence += line.strip()

                # Process the last sequence
                if sequence:
                    # Creating a temporary file for RNALFold output
                    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
                        temp_file_name = temp_file.name

                        # Call the RNALFold function, passing the sequence and the temp file name
                        # Replace this with actual call to RNALFold(sequence, temp_file_name)
                        RNA.Lfold(sequence, 150, temp_file)
                        # Read the output from the temporary file
                        temp_file.seek(0)
                        rnal_fold_output = temp_file.read()

                        min_value, min_index = extract_min_value_and_index(rnal_fold_output)
                        output_file.write(f">>{sequence_name}\nmax energy = {min_value} index = {min_index}\n\n")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # sequence = "GCCAATTGGTGCGGCAATTGATAATAACGAAAATGTCTTTTAATGATCTGGGTATAATGAGGAATTTTCCGAACGTTTTTACTTTATATATATATATACATGTAACATATATTCTATACGCTATAGAGAAAGGAAATTTT"  # Example long sequence
    # rnal_lowest_mfe_structure(sequence)
    process_RNA_sequences_v3('./input')

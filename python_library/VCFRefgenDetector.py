import gzip
import sys
import time
import pandas as pd
from VCF_header_refgen import *
from chromosomes_dict import *
from configobj import ConfigObj
from rich.console import Console
import argparse

final_results = []
config = ConfigObj("./config.config", interpolation="Template")
console = Console(highlight=False)

def gather_and_sum(lists):
    """It gathers and sums all the matches calculated in get_matches()"""
    cumulative_sums = {}
    for lst in lists:
        for key, value in lst:
            if key in cumulative_sums:
                cumulative_sums[key] += value
            else:
                cumulative_sums[key] = value
    console.print("Matches:", cumulative_sums)
    return cumulative_sums


def get_matches(snps, chr):
    """ Loads the pkl file depending on the chr that are present on the chunk. Compares the reference column from our
     input file and from the pkls, the number of matches to each version are returned to read_chunks()"""
    start = time.time()
    try:
        hg18 = pd.read_pickle(f"{config['pkl_paths']['pkl_hg18']}{chr}.pkl")
        GRCh37 = pd.read_pickle(f"{config['pkl_paths']['pkl_GRCh37']}{chr}.pkl")
        GRCh38 = pd.read_pickle(f"{config['pkl_paths']['pkl_GRCh38']}{chr}.pkl")
        t2t = pd.read_pickle(f"{config['pkl_paths']['pkl_T2T']}{chr}.pkl")
    except FileNotFoundError as e:
        console.print(f"FileNotFoundError: One or more pickle files for chromosome {chr} could not be found. Check the path "
              f"stated in the config file. Details: {e}")
        sys.exit(1)
    except (IOError, OSError) as e:
        console.print(f"IOError/OSError: An error occurred while reading the pickle files for chromosome {chr}. Details: {e}")
        sys.exit(1)


    console.print("Loading FP snps. Took:", time.time() - start, "s")

    start2 = time.time()
    genome_versions = [[hg18, "hg18"], [GRCh37, "GRCh37"], [GRCh38, "GRCh38"], [t2t, "T2T"]]
    matches = []

    for version in genome_versions:
        merged_df = version[0].set_index('position').join(snps.set_index('position'), how='inner')
        matches_count = (merged_df['nucleotide'] == merged_df[version[1]]).sum()
        matches.append([version[1], matches_count])

        # Uncomment if you want to get csv with the matches to each reference genome
        """matches_df = (merged_df['nucleotide'] == merged_df[version[1]]) 
        filtered_df = merged_df[matches_df]
        print(len(filtered_df))
        filtered_df.to_csv(f'matches_{version[1]}_{chr}.csv', mode='w')
        """

    # Display the DataFrame
    console.print("Getting matches. Took:", time.time() - start2, "s")

    return matches


def trimming_indels(content, ref):
    """If a row is longer than one position it is deleted, deleting this way any indels"""
    try:
        del_insertions = content[content.iloc[:, 2].str.len() == 1]  # delete insertions
        vcf_snps = del_insertions[del_insertions.iloc[:, 3].str.len() == 1]  # delete deletions
        snps = vcf_snps.iloc[:, [1, ref]]  # create df with pos and ref columns
        snps.iloc[:, 1] = snps.iloc[:, 1].apply(lambda x: x.upper())
        snps.columns = ['position', 'nucleotide']
        return snps
    except Exception:
        console.print("Reference column is empty, please check your input file. Stopping scan.")
        sys.exit(1)


def call_trimming(content, file_type, chr):
    """ The file must only contain SNPs. This function is necessary because the reference column has a different
    number in vcfs and in bim files """

    start = time.time()
    if file_type == "VCF":
        snps = trimming_indels(content, 2)  # content : chr pos ref alt
    elif file_type == "BIM":
        snps = trimming_indels(content, 3)
    console.print("Trimming indels. Took:", time.time() - start, "s")

    if len(snps) != 0:
        results = get_matches(snps, chr)
        final_results.append(results)
    else:
        console.print("There aren't FP SNPs in this chunk", style="bold")


def read_and_load(chunk, input_file, file_type):
    """To avoid loading a big pkl with information from all the chromosomes we first check which chr are there in
    the current chunk and then load only the necessary pkls
    """
    for chromosome, group_content in chunk.groupby(chunk.columns[0]):
        chromosome_str = str(chromosome)
        if chromosome_str in chromosome_map:
            chr_key = chromosome_map[chromosome_str]
            console.print("Variants being mapped from:", chr_key)
            call_trimming(group_content, file_type, chr_key)
        else:
            console.print(f"Chromosome {chromosome_str} not found in chromosome map. Skipping variants from {chromosome_str}.")


def read_chunks(complete_file, input_file, file_type, cols,  chunks, matches):
    """Loads the file in batches to avoid loading completely on memory. If there are enough matches to define the
    version, the loop will stop. If it has loaded more chunks than desired, the loop breaks too. For bim files,
    the inferred version must have at least 50% of the matches"""

    counter = 0
    results = {}
    console.print(
        f"Reading file in chunks of 100.000 variants - If file contains less than 100.000 variants the complete file "
        f"is used.")

    try:
        for chunk in pd.read_csv(complete_file, sep="\t", comment='#', header=None, chunksize=100000, usecols=cols):
            counter += 1 # TODO is this counter necessary?
            console.print("Chunk", counter)
            read_and_load(chunk, input_file, file_type)
            results = gather_and_sum(final_results)

            try:
                if max(results.values()) > matches:
                    break
                if counter > chunks:
                    break
            except ValueError:
                    console.print("0 FP SPNs in this chunk", style="bold")
            # Stop when the counter reaches 20
            if counter >= chunks:
                console.print(f"Stopping after processing {chunks} chunks.", style="bold")
                break
    finally:
        try:
            if max(results.values()) == 0:
                console.print("No SNPs found to infer the reference genome.", style="bold red")
            else:
                if max(results.values()) > sum(results.values())/2:
                    console.print("Inferred Reference genome:", max(results, key=results.get), style="bold dark_cyan")
                else:
                    console.print("Any of the versions have more than 50% of the total matches.  [bold] Reference "
                                  "genome version unknown. m", style="red")
        except ValueError as e:
            console.print("No SNPs found to infer the reference genome.", style="bold red")

def extract_columns(complete_file, input_file, file_type, chunks, matches):
    "Loads only the interesting columns"
    try:
        if file_type == "VCF":
            cols = [0, 1, 3, 4] # chr pos ref alt
            read_chunks(complete_file, input_file, file_type, cols, chunks, matches)
        elif file_type == "BIM":
            cols = [0, 3, 5, 4] # chr pos ref alt
            read_chunks(complete_file, input_file, file_type, cols, chunks, matches)
    except ValueError as e:
        console.print(f"ValueError: An error occurred while processing the columns for file {input_file}: {e}")
        sys.exit(1)
    except Exception as e:
        console.print(f"An unexpected error occurred while extracting columns from file {input_file}. Check your file "
              f"contains all mandatory columns. Error: {e}")
        sys.exit(1)


def extract_header(complete_file):
    "Extracts header and send it to match the refgenDetector database"
    header = []
    for line in complete_file:
        if line.startswith('#'):
            header.append(line.strip())
        else:
            break
    start_refgen_header(header)


def open_file(input_file, file_type, chunks, matches):
    """
    Checks if the file_type is correct. Then, opens the file if it endswith vcf/bim or vcf.gz/bim.gz.
    """
    console.print(f"Starting pre-processing for [[light_slate_blue]{input_file}[/light_slate_blue]]")
    formats = ("vcf", "bim")
    compressed_formats = ("vcf.gz", "bim.gz")
    dict_formats = {
        "vcf.gz": "VCF",
        "bim.gz": "BIM",
        "vcf": "VCF",
        "bim": "BIM"
    }

    for ext, expected_type in dict_formats.items():
        if input_file.endswith(ext):
            if file_type != expected_type:
                raise ValueError(f"File type should be {expected_type} for .{ext} files.")
            break

    try:
        if input_file.endswith(compressed_formats):
            with gzip.open(input_file, "rt") as complete_file:
                extract_header(complete_file)
                extract_columns(complete_file, input_file, file_type, chunks, matches)
        elif input_file.endswith(formats):
            with open(input_file, "rt") as complete_file:
                extract_header(complete_file)
                extract_columns(complete_file, input_file, file_type, chunks, matches)
        else:
            console.print("File format not supported. Please state a vcf or bim file, gz compressed or uncompressed")
    except FileNotFoundError:
        console.print(f"Error: The file {input_file} was not found in the stated path.")
        sys.exit(1)
    except OSError as e:
        console.print(f"Error: OS error occurred while handling the file {input_file}: {e}")
        sys.exit(1)
    except Exception as e:
        console.print(f"An unexpected error occurred while processing the file {input_file}: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(prog="INFERRING THE REFERENCE GENOME FROM A VARIANTS FILE")
    # MANDATORY ARGUMENTS
    parser.add_argument("-f", "--file", help="Path to the VCF or BIM file",
                        required=True)
    parser.add_argument("-t", "--type", choices=["VCF", "BIM"], help="Select the type of file stated in the --file "
                                                                     "argument",
                        required=True)

    # OPTIONAL ARGUMENTS
    parser.add_argument("-c", "--chunks", required=False, type=int, default=10000000000000000000000000,
                        help=("[OPTIONAL] By default the program will read the input file in chunks of 100.000 "
                              "variants. If you don't want to read the entire file select the maximum number of "
                              "chunks you " ## TODO delete this argument
                              "want to read. For example, -c 2 would read 200.000 variants, if the file has them."))
    parser.add_argument("-m", "--matches", required=False, type=int, default=5000,
                        help="[OPTIONAL] By the fault, when there are 5000 matches to a reference genome the reading "
                             "stops and results are print. You can modify the number of necessary matches with this "
                             "argument. ")  ## TODO modify the description : by default read the complete file ,
                                            # if used, set a number of matches to stop the run
    args = parser.parse_args()

    open_file(args.file, args.type, args.chunks, args.matches)


if __name__ == "__main__":
    main()



import gzip
import sys
import time
import pandas as pd
from dns.inet import inet_pton
from VCF_header_refgen import *
from chromosomes_dict import *
from configobj import ConfigObj
from rich.console import Console
import json

final_results = []
config = ConfigObj("./config.config", interpolation="Template")
console = Console(highlight=False)

# Initialize an empty dictionary for the JSON structure
result_dict = {
    "file_name": None,
    "VCFRefgenDetector": {
        "info_from_header": {
            "inferred_reference_genome": None,
            "gVCF": None,
            "n_of_samples": None
        },
        "info_from_variants": {
            "inferred_reference_genome": None,
            "gVCF": None
        },
        "date": None,
        "version": None
    }
}

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

    return matches


def trimming_indels(content, ref):
    """If a row is longer than one position it is deleted, deleting this way any indels"""
    #try:
    del_insertions = content[(content.iloc[:, 2].str.len() == 1) | (content.iloc[:, 2] == '<NON_REF>')]
    vcf_snps = del_insertions[
        (del_insertions.iloc[:, 3].str.len() == 1) | (del_insertions.iloc[:, 3] == '<NON_REF>')]
    snps = vcf_snps.iloc[:, [1, ref]]  # create df with pos and ref columns
    snps.iloc[:, 1] = snps.iloc[:, 1].apply(lambda x: x.upper())
    snps.columns = ['position', 'nucleotide']
    if content.iloc[:, 3].astype(str).str.contains('<NON_REF>').any():
        print("\t gVCF by ALT column")
        result_dict["VCFRefgenDetector"]["info_from_variants"]["gVCF"] = "True"
    else:
        result_dict["VCFRefgenDetector"]["info_from_variants"]["gVCF"] = "False"
    return snps
"""    except Exception:
        console.print("Reference column is empty, please check your input file. Stopping scan.")
        sys.exit(1)"""


def call_trimming(content, chr):
    """ The file must only contain SNPs. This function is necessary because the reference column has a different
    number in vcfs and in bim files """

    snps = trimming_indels(content, 2)  # content : chr pos ref alt

    if len(snps) != 0:
        results = get_matches(snps, chr)
        final_results.append(results)
    else:
        console.print("There aren't FP SNPs in this chunk", style="bold red")


def read_and_load(chunk):
    """To avoid loading a big pkl with information from all the chromosomes we first check which chr are there in
    the current chunk and then load only the necessary pkls
    """
    console.print("[bold]Variants info: [/bold]")
    for chromosome, group_content in chunk.groupby(chunk.columns[0]):
        chromosome_str = str(chromosome)
        if chromosome_str in chromosome_map:
            chr_key = chromosome_map[chromosome_str]
            console.print("Variants being mapped from:", chr_key)
            call_trimming(group_content, chr_key)
        else:
            console.print(f"Chromosome {chromosome_str} not found in chromosome map. Skipping variants from {chromosome_str}.")


def read_chunks(complete_file, cols):
    """Loads the file in batches to avoid loading completely on memory. If there are enough matches to define the
    version, the loop will stop. If it has loaded more chunks than desired, the loop breaks too. For bim files,
    the inferred version must have at least 50% of the matches"""
    counter = 0
    results = {}

    try:
        for chunk in pd.read_csv(complete_file, sep="\t", comment='#', header=None, chunksize=100000, usecols=cols):
            counter += 1 # TODO is this counter necessary?
            read_and_load(chunk)
            results = gather_and_sum(final_results)

            try:
                if counter > int(config['chunks']['n_chunks']):
                    break
            except ValueError:
                    console.print("0 FP SPNs in this chunk", style="bold")
            # Stop when the counter reaches 20
            if counter >= int(config['chunks']['n_chunks']):
                break
    finally:
        try:
            if max(results.values()) == 0:
                console.print("No SNPs found to infer the reference genome.", style="bold red")
                result_dict["VCFRefgenDetector"]["info_from_variants"]["inferred_reference_genome"] = "No SNPs found to infer the reference genome."
            else:
                if max(results.values()) > sum(results.values())/2:
                    console.print("\t Inferred Reference genome:", max(results, key=results.get))
                    result_dict["VCFRefgenDetector"]["info_from_variants"][
                        "inferred_reference_genome"] = max(results, key=results.get)

                else:
                    console.print("Some of the versions have more than 50% of the total matches.  [bold] Reference "
                                  "genome version unknown. m", style="red")
                    result_dict["VCFRefgenDetector"]["info_from_variants"][
                        "inferred_reference_genome"] = "Reference genome unknown"


            result_dict["VCFRefgenDetector"]["date"] = time.strftime("%Y-%m-%d")
            result_dict["VCFRefgenDetector"]["version"] = config["version"]["version"]

        except ValueError as e:
            console.print("No SNPs found to infer the reference genome.", style="bold red")

def extract_columns(complete_file):
    "Loads only the interesting columns"
    cols = [0, 1, 3, 4] # chr pos ref alt
    read_chunks(complete_file, cols)

def get_n_samples(header):

    mandatory_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] #fixed fields for variant information
    columns = header[-1].split("\t")

    for column in mandatory_columns:
        columns.remove(column)

    result_dict["VCFRefgenDetector"]["info_from_header"]["n_of_samples"] = len(columns)


def extract_header(complete_file):
    "Extracts header and send it to match the refgenDetector database"
    header = []
    for line in complete_file:
        if line.startswith('#'):
            header.append(line.strip())
        else:
            break

    get_n_samples(header)

    gVCF, ref_gen = start_refgen_header(header)

    result_dict["VCFRefgenDetector"]["info_from_header"]["inferred_reference_genome"] = ref_gen

    if gVCF:
        result_dict["VCFRefgenDetector"]["info_from_header"]["gVCF"] = "True"
    else:
        result_dict["VCFRefgenDetector"]["info_from_header"]["gVCF"] = "False"

def open_file(input_file):

    result_dict["file_name"] = input_file

    console.print(f"[bold]File: [/bold] {input_file}")
    formats = ("vcf", "bim")
    compressed_formats = ("vcf.gz", "bim.gz")

    if input_file.endswith(compressed_formats):
        with gzip.open(input_file, "rt") as complete_file:
            extract_header(complete_file)
            extract_columns(complete_file)

    elif input_file.endswith(formats):
        with open(input_file, "rt") as complete_file:
            extract_header(complete_file)
            extract_columns(complete_file)

    # Once all the results are gathered, save the dictionary to a JSON file
    with open(config["json"]["path_json"], "a") as json_file:
        json.dump(result_dict, json_file, indent=4)

open_file(sys.argv[1])





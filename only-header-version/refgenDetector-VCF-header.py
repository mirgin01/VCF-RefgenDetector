"""
refgenDetector sub-version
10/12/2025
Mireia Marin i Ginestar
This refgenDetector-VCF-Header infers the reference genome ONLY from the header of the VCF file.
*Important Note* the contig information is NOT mandatory for VCFs. The reference genome not being inferred doesn't mean there's something missing or incorrect in the VCF file. 
"""

import gzip
import sys
from aligment_files import comparison
from chromosomes_dict import *
from rich.console import Console
from refgenDetector.reference_genome_dictionaries import *

console = Console(highlight=False)

def intersection_targetfile_referencerepo(dict_SN_LN, reference_genome):
    """
    Find the matches between the target file and the repository of unique contigs per reference genome.
    Returns the actual matches (ln info) instead of just their count.
    Args:
         dict_SN_LN (dict) : dictionary with the contig (SN: key, LN: value) info from the target file
         reference_genome (dict entry): one of the versions from major_releases

    Returns:
        matches (set) : list of lengths matching to the version currently being read
        reference_genome["build"] (str): build of the version currently being read
        reference_genome["species"] (str): species of the version currently being read
    """
    matches = set(dict_SN_LN.values()).intersection(reference_genome["ref_gen"].values())
    return matches, reference_genome["build"], reference_genome["species"]

def check_if_decoy(matches_info, target_file): 
    """
    Checks if there's inconsistency of the versions in the target file or if the multiple matches are random
    Args:
         matches_info (list): list of tuples. Each tuple have 3 positions: lengths from the contigs matching,
         version where the contigs match, species from the version.
         target_file (str): path of the target file

    Returns:
        If the matches to the secondary version are at least as long as the shortest chromosome of the version with
        more matches then a message raising the incosistency is printed.
        If the matches to the secondary version are shorter than the shortest chromosome then it assumes it's a decoy
        contig matching another version randomly and returns:
        incosistency = False so the code can continue in comparison() and give the results based on the version with
        most matches.
    """
    incosistency_found = False
    filtered_entries = [entry for entry in matches_info if entry[0]]  # get the multiple matches
    match = max(filtered_entries, key=lambda ref_gen_w_macthes: len(ref_gen_w_macthes[0]))  # version with most
    # matches with LN values
    inconsistency_matches = []
    for version in filtered_entries:
        if version != match:
            for ln in version[0]:
                if int(ln) > min_values[match[1]]: # checks if the ln matching is more or less chr length
                    ref_dict = globals().get(version[1])
                    inconsistency_matches.extend([key for key, value in ref_dict.items() if value == ln])
                    incosistency_found = True
            if incosistency_found == True:
                console.print(f"[bold]File:[/bold] {target_file} \n[bold][red]Error:[/bold] Inconsistency found "
                              f"- file contains contigs from different genome versions[/red]")
                console.print(f"[red]Contigs {inconsistency_matches} belong to {version[1]}, but the rest belongs to"
                              f" {match[1]}[/red].")
    return incosistency_found

def comparison(dict_SN_LN, target_file):
    """
    First, it defines the major release to which the header belongs to. Then, checks if a flavor can be inferred.
    Args:
         dict_SN_LN (dict): dictionary with the contig (SN: key, LN: value) info from the target file
         target_file (str): path of the target file

    Returns:
        Prints the file path being analyzed,the species and the reference genome version inferred.
        It raises an error if:
            - The contigs in the target file are not in the database (a species or ref gen version not included in the tool)
            - There are contigs belonging to more than one release/species. This will be printed if the match between
            species is as long as the shortest chromosome from the version with the most matches. If the match is
            shorter it assumes it's a random match e.g a decoy contig that randomly matches the length of
            another species/version.
    """

    matches_info = [intersection_targetfile_referencerepo(dict_SN_LN, major_releases[ref]) for ref in major_releases]
    matches_with_counts = [(len(matches), build, species) for matches, build, species in matches_info]
    max_match = max(matches_with_counts, key=lambda ref_gen_w_macthes: ref_gen_w_macthes[0]) # Find the major release
    # with the maximum matches
    incosistency = False 

    # check all the matches belong to the same release version
    multiple_matches = []
    for match in matches_with_counts:
        if match[0] != 0:
            multiple_matches.append(match)

    if len(multiple_matches) > 1 :
        if multiple_matches[0][1] != "hg17" and multiple_matches[1][1] != "hg18": # these versions share contig lengths
            if multiple_matches[0][1] != "rhemac3" and multiple_matches[1][1] != "rhemac8":
                incosistency = check_if_decoy(matches_info, target_file)

    if incosistency == False:

        if max_match[0] == 0:
            console.print(f"[bold]File:[/bold] {target_file} \n[bold][red]Reference genome can't be inferred[/bold] - "
                          "The contigs in the file are not found in refgenDetector database[red]")

        elif max_match[1] == "GRCh37": #check for GRCh37 flavors

            matches_flavors = [
                intersection_targetfile_referencerepo(dict_SN_LN, flavors_GRCh37[ref])
                for ref in flavors_GRCh37
            ]
            
            match_flavors = max(matches_flavors, key=lambda x: x[0])
            if match_flavors: #if some flavor was defined it prints it
                console.print(f"[bold]File:[/bold] {target_file} \n[bold]Species detected:[/bold] {match_flavors[2]} "
                f"[bold]\nReference genome version:[/bold] {match_flavors[1]}")
            else: #if there wasnt any flavor inferred, the major release it printed
                console.print(f"[bold]File:[/bold] {target_file} \n[bold]Species detected:[/bold] Homo sapiens \n["
                              f"bold]Reference genome version:[/bold] GRCh37")

        elif max_match[1] == "GRCh38": #checks for GRCh38 flavors

            if any("HLA-" in key for key in dict_SN_LN.keys()):
                #first checks if the contigs contain in their names HLA-
                console.print(f"[bold]File:[/bold] {target_file} \n[bold]Species detected:[/bold] Homo sapiens \n[bold]"
                              f"Reference genome version:[/bold] hs38DH_extra")
            elif set(dict_SN_LN.values()).intersection(verily_difGRCh38.values()):#checks if the Verily's unique
                # lengths are present
                console.print(f"[bold]File:[/bold] {target_file} \n[bold]Species detected:[/bold] Homo sapiens \n[bold]"
                              f"Reference genome version:[/bold] GRCh38_no_alt_plus_hs38d1")
            else: # if no GRCh38 flavor is inferred, the major release is printed
                console.print(f"[bold]File:[/bold] {target_file} \n[bold]Species detected:[/bold] Homo sapiens \n["
                              f"bold]Reference genome version:[/bold] GRCh38")
        else: # print the major releases with no considered flavors.
            console.print(f"[bold]File:[/bold] {target_file} \n[bold]Species detected:[/bold] {match[2]} "
                  f"\n[bold]Reference genome version:[/bold] {match[1]}")

def start_refgen_header(header):
    contig_list = [line for line in header if '##contig' in line and 'length' in line]
    if len(contig_list)!=0:
        contig_list2= [i.split(",") for i in contig_list]
        dict_contigs = {}

        for line in contig_list2:
            contig_id = None
            contig_length = None

            for part in line:
                if part.startswith('##contig=<ID='):
                    contig_id = part.replace('##contig=<ID=', '')
                elif part.startswith('length='):
                    contig_length = int(part.replace('length=', '').replace('>', ''))

            if contig_id is not None and contig_length is not None:
                dict_contigs[contig_id] = contig_length
        comparison(dict_contigs, sys.argv[1])  # run the next f

    else:
        console.print("[dark_orange]Contig information not in the header[/dark_orange] - [bold dark_orange]The reference genome can't be "
                      "inferred from the header information [/bold dark_orange]")
        
def extract_header(complete_file):
    "Extracts header and send it to match the refgenDetector database"
    header = []
    for line in complete_file:
        if line.startswith('#'):
            header.append(line.strip())
        else:
            break
    start_refgen_header(header)


def main():
    try:  # tries to open the file with utf-8 encoding
        with open(sys.argv[1], "r") as header_txt:
            extract_header(header_txt)
    except UnicodeError:  # tries to open the file with iso-8859-1 encoding
        with open(sys.argv[1], "r", encoding="iso-8859-1") as header_txt:  # tries to open the file with utf-8
            # encoding
            extract_header(header_txt)

if __name__ == "__main__":
    main()
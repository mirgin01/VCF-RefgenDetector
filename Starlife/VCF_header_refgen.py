from reference_genome_dictionaries import GRCh38, GRCh37, hs37d5, hg16, hg17, hg18, hg19, b37, verily_difGRCh38, T2T
import sys
from rich.console import Console

console = Console(highlight=False)


def intersection_targetfile_referencerepo(dict_SN_LN, reference_genome):
    """
    Find the matches between the target file and the repository of unique contigs per reference genome
    """
    return set(dict_SN_LN.values()).intersection(reference_genome.values())

def comparison(dict_SN_LN):

    """
    First, it defines the major release to which the header belongs to. Then, checks if there's any match with the flavors.
    """
    major_releases = {"hg16": hg16,"hg17": hg17,"hg18": hg18,"GRCh37": GRCh37,"GRCh38": GRCh38,"T2T": T2T} #major
    # release that can be inferred
    flavors_GRCh37 = {"hs37d5": hs37d5,"b37": b37,"hg19": hg19} #GRCh37 flavors that can be inferred
    major_release_list = [major_releases[ref] for ref in major_releases if
                          intersection_targetfile_referencerepo(dict_SN_LN, major_releases[ref])] #gets the major
    # release to which the header belongs to
    if len(major_release_list) == 0: # if there wasnt any major release inferred print:
        console.print( "\t Reference genome can't be inferred")
        ref_gen = "Reference genome can't be inferred"
    else:
        if major_release_list[0] == major_releases["GRCh37"]: #check for GRCh37 flavors
            match_flavor = next(
                (flavors_GRCh37[flav] for flav in flavors_GRCh37 if intersection_targetfile_referencerepo(dict_SN_LN,
                                                                                             flavors_GRCh37[flav])),
                None) #infers the flavor of GRCh37 to which the header belongs to
            if match_flavor: #if some flavor was defined it prints it
                console.print(f"\t Reference genome: "
                      f":{[k for k, v in flavors_GRCh37.items() if v == match_flavor][0]}")
                ref_gen = f"{[k for k, v in flavors_GRCh37.items() if v == match_flavor][0]}"
            else: #if there wasnt any flavor inferred, the major release it printed
                console.print("\t Reference genome: GRCh37")
                ref_gen = "GRCh37"
        elif major_release_list[0] == major_releases["GRCh38"]: #checks for GRCh38 flavors
            if any("HLA-" in key for key in dict_SN_LN.keys()): #first checks if the contigs contain in their names HLA-
                console.print("\t Reference genome: hs38DH_extra") #if so,
                # the reference genome
                # used was hs38DH_extra
                ref_gen = "hs38DH_extra"
            elif intersection_targetfile_referencerepo(dict_SN_LN, verily_difGRCh38):#checks if the Verily's unique
                # lengths are present
                console.print("\t Reference genome:Verily's GRCh38")
                ref_gen = "Verily's GRCh38"
            else: # if no GRCh38 flavor is inferred, the major release is printed
                console.print("\t Reference genome: GRCh38")
                ref_gen = "GRCh38"
        else: #print the major releases with no considered flavors.
            console.print(f"\t Reference genome:"
                  f"{[k for k, v in major_releases.items() if v == major_release_list[0]][0]}")
            ref_gen = f"{[k for k, v in major_releases.items() if v == major_release_list[0]][0]}"

    return ref_gen


def start_refgen_header(header):
    console.print("[bold]Header info: [/bold]")

    gVCF= [line for line in header if '##ALT=<ID=NON_REF' in line]
    if gVCF:
        print("\tgVCF according to header")

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

        ref_gen = comparison(dict_contigs)

    else:
        console.print("\tContig information not in the header")
        ref_gen = "Contig information not in the header"

    return gVCF , ref_gen





# VCF-RefgenDetector

## How to install 

1. Clone the repository
2. Get the pkl folder

## How to use
```
[bio-box: /bio-scratch/mireia/VCF-RefgenDetector/VCF-RefgenDetector] # python VCFRefgenDetector.py -h
usage: INFERRING THE REFERENCE GENOME FROM A VARIANTS FILE [-h] -f FILE -t {VCF,BIM} [-c CHUNKS] [-m MATCHES]

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Path to the VCF or BIM file
  -t {VCF,BIM}, --type {VCF,BIM}
                        Select the type of file stated in the --file argument
  -c CHUNKS, --chunks CHUNKS
                        [OPTIONAL] By default the program will read the input file in chunks of 100.000 variants. If you don't want to read the entire file select the maximum number of chunks you want
                        to read. For example, -c 2 would read 200.000 variants, if the file has them.
  -m MATCHES, --matches MATCHES
                        [OPTIONAL] By the fault, when there are 5000 matches to a reference genome the reading stops and results are print. You can modify the number of necessary matches with this
                        argument.
```

## Example

```
[bio-box: /bio-scratch/mireia/VCF-RefgenDetector/VCF-RefgenDetector] # python VCFRefgenDetector.py -f ../1000VCFs/10.vcf.gz -t VCF
Starting pre-processing for [../1000VCFs/10.vcf.gz]
The reference genome can't be inferred from the header information 
Reading file in chunks of 100.000 variants - If file contains less than 100.000 variants the complete file is used.
Chunk 1
Variants being mapped from: chr10
Trimming indels. Took: 0.08820486068725586 s
Loading FP snps. Took: 2.690643787384033 s
Getting matches. Took: 8.830072164535522 s
Matches:
{'hg18': 1, 'GRCh37': 9230, 'GRCh38': 2, 'T2T': 2}
Inferred Reference genome: GRCh37
```

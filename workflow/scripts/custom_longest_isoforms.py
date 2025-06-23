## this script is a modification of the primary_transcript.py file included
## in David Emms's Orthofinder software

import os
import re
import sys
import argparse
from collections import Counter, defaultdict

# Use the 'all' version rather than ab initio

def CheckFile(fn):
    """
    Checks for:
    - Duplicated accession lines
    """
    accs = set()
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith(">"):
                a = l.rstrip()[1:]
                if a in accs:
                    print("\nERROR: duplicated sequence accession:\n%s" % a)
                    print("\nPlease correct this and then rerun the script.\n")
                    return False
                accs.add(a)
    return True

def ScanTags(fn):
    """
    For ensembl genomes, look for tag:id and count repeated ids
    :param fn:
    :return:
    """
    tags = set()
    tokens = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            tokens.append([t.split(":", 1) for t in line.rstrip().split() if ":" in t])
            tags.update([t[0] for t in tokens[-1]])
    for this_tag in tags:
        print(this_tag)
        # print(tokens[-1])
        c = Counter([idd for acc in tokens for t, idd in acc if t == this_tag])
        print(c.most_common(5))
        print("")

def ScanTags_NCBI(fn):
    genes = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            genes.append(line[1:].split(".", 1)[0])
    print("%d sequences, %d genes" % (len(genes), len(set(genes))))

def ScanTags_with_fn(fn, gene_name_fn):
    genes = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            genes.append(gene_name_fn(line))
    print("%d sequences, %d genes" % (len(genes), len(set(genes))))
    # print(genes[0])
    # print(sorted(genes)[:10])

def GetGeneName_Ensembl(acc_line):
    tokens = [(t.split("=") if "=" in t else t.split(":"))[1] for t in acc_line.rstrip().split() if ("gene:" in t or "gene=" in t or "locus:" in t or "locus=" in t)]
    if len(tokens) != 1: return None
    return tokens[0]

def GetGeneName_Custom(acc_line, keywords, prefix="", suffix=""):
    tokens = [
        (t.split("=") if "=" in t else t.split(":"))[1]
        for t in acc_line.rstrip().split()
        if any(f"{keyword}:" in t or f"{keyword}=" in t for keyword in keywords)
    ]
    #print(tokens)
    if len(tokens) != 1:
        return None
    return f"{prefix}{tokens[0]}{suffix}"

def GetGeneName_ByKey(acc_line, key="gene"):
    """
    Extracts the value after 'key=' in the header line.
    Returns None if not found.
    """
    pattern = rf'{re.escape(key)}=([^; ]+)'
    match = re.search(pattern, acc_line)
    if match:
        return match.group(1)
    return None

def IsNCBI(fn):
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith(">"):
                l = l.rstrip()
                if l.startswith(">NP_") and l.endswith("]"): return True
                elif l.startswith(">XP_") and l.endswith("]"): return True
                elif l.startswith(">YP_") and l.endswith("]"): return True
                elif l.startswith(">WP_") and l.endswith("]"): return True
                return False
    return False

def GetGeneName_NCBI(acc_line):
    acc_line = acc_line[1:]
    original = acc_line
    # look for "isoform X[:d]+" or "isoform [:d]+"
    acc_line = re.sub("isoform [0-9, A-Z]+ ", "", acc_line)
    acc_line = re.sub("isoform X[0-9, A-Z]+ ", "", acc_line)
    # This last step is nasty. These are the same gene:
    # >XP_024356342.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
    # >XP_024356343.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
    # as the name is the same and they same 'isoform ...'
    # Whereas these are not the same gene, even though the names are identical
    # because they don't say isoform:
    # >XP_024390255.1 40S ribosomal protein S12-like [Physcomitrella patens]
    # >XP_024399722.1 40S ribosomal protein S12-like [Physcomitrella patens]
    #
    # To deal with that, we remove the ID (e.g. XP_024356342.1) if it says 'isoform'
    # so that the lines are identical, but not when it doesn't say 'isoform'
    # so that the lines are different. If I were writting the script from scratch
    # for NCBI files I'd do it a different way, but this is a way to handle it so 
    # that it works with the existing logic in the file.
    if original != acc_line:
        acc_line = acc_line.split(None, 1)[-1]
    return acc_line

def CreatePrimaryTranscriptsFile(fn, dout, gene_name_fn, q_use_original_accession_line, remove_pattern=None, prefix="", suffix="", **kwargs):
    # Get genes and lengths
    max_gene_lens = defaultdict(int)
    acc_to_use = defaultdict(str)
    nAcc = 0
    nGeneUnidentified = 0

    with open(fn, 'r') as infile:
        lines = [l.rstrip() for l in infile]
    N = len(lines) - 1

    # First pass: Identify the longest sequence for each gene
    iLine = -1
    while iLine < N:
        iLine += 1
        line = lines[iLine]
        if not line.startswith(">"):
            continue
        nAcc += 1
        iLineAcc = iLine
        gene = gene_name_fn(line, **kwargs) if kwargs else gene_name_fn(line)
        if gene is None:
            # Fallback: Use the full header (minus ">") if gene_name_fn fails
            gene = line[1:]
            if remove_pattern:
                gene = gene.replace(remove_pattern, "")
            nGeneUnidentified += 1
        else:
            # Apply remove_pattern to the extracted gene name
            if remove_pattern:
                gene = gene.replace(remove_pattern, "")

        # Calculate sequence length
        l = 0
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if line.startswith(">"):
                iLine -= 1
                break
            l += len(line.rstrip())

        # Update the longest sequence for the gene
        if l > max_gene_lens[gene]:
            max_gene_lens[gene] = l
            acc_to_use[gene] = iLineAcc

    print("Found %d accessions, %d genes, %d unidentified transcripts" % (nAcc, len(max_gene_lens), nGeneUnidentified))

    # Second pass: Write the longest sequences to the output file
    nGenesWriten = 0
    with open(dout, 'w') as outfile:
        iLine = -1
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if not line.startswith(">"):
                continue
            gene = gene_name_fn(line, **kwargs) if kwargs else gene_name_fn(line)
            if gene is None:
                # Fallback: Use the full header (minus ">") if gene_name_fn fails
                gene = line[1:]
                if remove_pattern:
                    gene = gene.replace(remove_pattern, "")
            else:
                # Apply remove_pattern to the extracted gene name
                if remove_pattern:
                    gene = gene.replace(remove_pattern, "")

            # Skip sequences that are not the longest for their gene
            if gene not in acc_to_use or iLine != acc_to_use[gene]:
                continue

            # Add prefix and suffix to the gene name
            gene = f"{prefix}{gene}{suffix}"

            # Write the header
            if q_use_original_accession_line:
                # Use the original header as-is
                acc_line_out = line + "\n"
            elif gene_name_fn == GetGeneName_NCBI:
                # For NCBI, only use the ID part of the header (everything before the first space)
                acc_line_out = line.split(" ", 1)[0] + "\n"
            else:
                # Default behavior: Use the extracted gene name with prefix and suffix
                acc_line_out = f">{gene}\n"
            outfile.write(acc_line_out)

            # Write the sequence
            nGenesWriten += 1
            while iLine < N:
                iLine += 1
                line = lines[iLine]
                if line.startswith(">"):
                    iLine -= 1
                    break
                outfile.write(line + "\n")

    print("Wrote %d genes" % nGenesWriten)
    if nGenesWriten != len(max_gene_lens):
        print("ERROR: Mismatch in the number of genes written.")
        raise Exception
    print("Output written to:", dout)


def last_dot(text):
    return text[1:].rstrip().rsplit(".", 1)[0]


def space(text):
    return text[1:].rstrip().split(None, 1)[0]


function_dict = {"last_dot": last_dot, "space": space}

def main():
    parser = argparse.ArgumentParser(description="Extract longest isoforms from a FASTA file.")
    parser.add_argument("input_file", help="Path to the input FASTA file.")
    parser.add_argument("output_file", help="Path to the output FASTA file.")
    parser.add_argument("--gene_name_function", choices=["ensembl", "ncbi", "customized", "by_key"], default="ensembl",
                        help="Function to use for extracting gene names. Default is 'ensembl'.")
    parser.add_argument("--key", default="gene",
                        help="Key to search for in the header when using 'by_key' function (default: 'gene').")
    parser.add_argument("--keywords", nargs="*", default=[], help="name of the gene id to use for filtering for the custom gene name function.")
    parser.add_argument("--prefix", default="", help="Prefix to add to gene names (used with 'custom').")
    parser.add_argument("--suffix", default="", help="Suffix to add to gene names (used with 'custom').")
    parser.add_argument("--remove_pattern", default="", help="pattern to remove from gene ids).")
    parser.add_argument("--use_original_accession", action="store_true", default=False,
                        help="Flag to use the original accession line in the output. Default is False.")

    args = parser.parse_args()

    if not CheckFile(args.input_file):
        return

    gene_name_fn = None
    kwargs = {}

    if args.gene_name_function == "ensembl":
        gene_name_fn = GetGeneName_Ensembl
    elif args.gene_name_function == "ncbi":
        gene_name_fn = GetGeneName_NCBI
    elif args.gene_name_function == "custom":
        gene_name_fn = GetGeneName_Custom
        kwargs = {"keywords": args.keywords, "prefix": args.prefix, "suffix": args.suffix}
    elif args.gene_name_function == "by_key":
        gene_name_fn = GetGeneName_ByKey
        kwargs = {"key": args.key}

    CreatePrimaryTranscriptsFile(
        args.input_file,
        args.output_file,
        gene_name_fn,
        args.use_original_accession,
        args.remove_pattern,
        prefix=args.prefix,
        suffix=args.suffix,
        **kwargs
    )


if __name__ == "__main__":
    main()

Mongrail.2.0 Documentation

# Installation


## Compiling the program (Linux)

Download the file [mongrail2.v1.0.0.tar.gz] and unzip the archive. The source code of the program is written in C. The pkg-config and glib-2.0 libraries must be installed prior to compiling. You'll need a C compiler installed. The program requires bcftools and tabix to be installed. The steps for compiling are as follows:

```
tar -xvzf *.tar.gz
cd mongrail2*
make all
```

This will create the executable files `mongrail2` and `gendiplo` in the current directory and example data files (in subfolder example).

# Running the program

## Getting Mongrail up and running on Mac OS X or Unix

The subfolder example contains a puma dataset (34 scaffolds with 10 markers each). To run the puma dataset type:

```
./run_mongrail.sh -A example/Texas.vcf.gz -B example/Cypress.vcf.gz -i example/hybrids.vcf.gz -r 1.9
```

# Input files

The mongrail program requires three input files:

1. *phasedA* file - a Variant Call Format (VCF) file containing the phased data for population A.

2. *phasedB* file - a Variant Call Format (VCF) file containing the phased data for population B.

3. *hybrid* file - a Variant Call Format (VCF) file containing the multi-locus genotype data for hybrids.

Currently the program can handle a maximum of 10 markers per chromosome. If the VCF file contains more than 10 markers per chromosome, the program will choose the first 10 markers. 

# Required options and arguments

The mongrail2 program requires the following three options and arguments:

| Option | Value    | Description                    |
|--------|----------|--------------------------------|
| \-A    | filename | Population A VCF file (phased) |
| \-B    | filename | Population B VCF file (phased) |
| \-i    | filename | Hybrids VCF file               |

# Command line options

The mongrail2 program needs you to specify the recombination rate (in cM/Mb) to be considered for each chromosomes. If the recombination rate is same for all chromosomes, use option `-r` followed by the value. If you want to use different rates for different chromosomes, use option `-R` followed by the filename that contains the chromosome name and corresponding recombination rates. For example, each line of the input file should have the following format

`chromID rec_rate`

where `chromID` is the name of the chromosome and `rec_rate` is the corresponding recombination rate (in cM/Mb).

# Interpreting the output files

*output.txt* is a text file that contains the estimated posterior probabilities that each individual (input from the *hybrid* file) belongs to each of the six different genealogical classes. The classes are ordered as follows: a, b, c, d, e, f.

*post_prob.txt* is a text file that gives a detailed description of the estimated posterior probabilities for each individual under every chromosome. 

*output.txt* file is nothing but a summary of *post_prob.txt*. It contains the final posterior probabilities calculated over all chromosomes.

# MONGRAIL 2.0

Bayesian inference of hybrid genealogy from population genetic data.

MONGRAIL classifies individuals as purebred or hybrid based on multilocus genotype data from two reference populations. It computes log-likelihoods for six genealogical models representing different hybrid categories.

## Installation

### Requirements

- GCC compiler
- GNU Make
- Standard C math library
- htslib (for VCF support)

### Build

```bash
git clone https://github.com/mongrail/mongrail2.git
cd mongrail2/brucesrc
make
```

This produces the `mongrail2` executable.

### Build options

The default build uses `-O3 -march=native -ffast-math` for optimal performance. For debugging, edit the Makefile:

```makefile
CFLAGS = -Wall -g
```

## Usage

```
mongrail2 [options] popA popB hybrids
```

### Arguments

| Argument | Description |
|----------|-------------|
| `popA` | Genotype file for reference population A |
| `popB` | Genotype file for reference population B |
| `hybrids` | Genotype file for putative hybrid individuals |

### Options

| Option | Description |
|--------|-------------|
| `-c` | Read VCF format (requires phased, biallelic SNPs; max 16 per chromosome) |
| `-r RATE` | Recombination rate per base pair (default: 1e-8) |
| `-v` | Verbose output (show input file summary) |
| `-h` | Show help message |
| `-V` | Show version |

### Example

```bash
# Using .GT format
./mongrail2 popA.GT popB.GT hybrids.GT

# Using VCF format
./mongrail2 -c popA.vcf popB.vcf hybrids.vcf

# With options
./mongrail2 -c -r 1e-7 -v popA.vcf popB.vcf hybrids.vcf
```

## Input File Format

Input files use the `.GT` format with tab-separated fields:

```
ChromName  Position  Ind1  Ind2  Ind3  ...
```

Each individual's genotype is encoded as `allele1|allele2` or `allele1/allele2`:
- `|` indicates phased data
- `/` indicates unphased data
- Alleles are `0` or `1`

### Example input file

```
Chr1    1000    0|0    0|1    1|1
Chr1    2000    0|1    1|1    0|0
Chr1    3500    0|0    0|1    0|1
Chr2    500     1|0    0|0    1|1
Chr2    1200    0|1    0|1    0|0
```

### Constraints (.GT format)

- Maximum 16 SNPs per chromosome/scaffold
- Maximum 1000 chromosomes/scaffolds
- Maximum 1000 individuals per population
- All three input files must have the same chromosomes and SNP positions

## VCF Input Format

When using the `-c` option, MONGRAIL reads standard VCF files.

### Requirements

- **Phased genotypes**: All genotypes must use `|` separator (e.g., `0|1`)
- **Biallelic SNPs only**: No indels or multiallelic sites
- **Maximum 16 SNPs per chromosome**: Hard limit for computational tractability
- **Matching positions**: All three VCF files must have identical chromosomes and positions

### Example VCF

```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS   ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  IND1  IND2  IND3
chr1    1000  .   A    G    .     PASS    .     GT      0|0   0|1   1|1
chr1    2000  .   C    T    .     PASS    .     GT      0|1   1|1   0|0
chr1    3500  .   G    A    .     PASS    .     GT      0|0   0|1   0|1
```

### Error messages

| Error | Cause |
|-------|-------|
| `variant at X:Y is not a biallelic SNP` | Indel or multiallelic variant found |
| `unphased or missing genotype at X:Y` | Genotype uses `/` instead of `|`, or is missing |
| `too many SNPs on chromosome X` | More than 12 SNPs on a single chromosome |

## Output Format

Output is written to stdout in a format suitable for downstream processing:

```
# mongrail 2.0
# recomb_rate=1.00e-08  chroms=2  n_A=4  n_B=4  n_hyb=3
#
# indiv      Ma          Mb          Mc          Md          Me          Mf
      0     -0.8825     -2.2860     -9.0206    -12.9229    -10.3721     -3.6376
      1    -12.9229    -10.2429     -8.8892     -0.6351     -2.0255     -3.3791
      2     -7.4921     -2.1017     -0.7920     -7.4921     -2.1017     -2.1018
```

- Header lines begin with `#`
- Each row contains log-likelihoods for one individual across all six models
- The model with the highest (least negative) log-likelihood is the best fit

### Filtering output

```bash
# Data only (no header)
./mongrail2 popA.GT popB.GT hybrids.GT | grep -v '^#'

# Pipe to file
./mongrail2 popA.GT popB.GT hybrids.GT > results.txt
```

## Genealogical Models

MONGRAIL evaluates six genealogical models for each individual:

| Model | Code | Description | Parental Origin |
|-------|------|-------------|-----------------|
| **Ma** | d | Pure population A | Both haplotypes from pop A |
| **Mb** | b | Backcross to A | One haplotype from A, one recombinant (F1 parent) |
| **Mc** | c | F1 hybrid | One haplotype from A, one from B |
| **Md** | a | Pure population B | Both haplotypes from pop B |
| **Me** | e | Backcross to B | One haplotype from B, one recombinant (F1 parent) |
| **Mf** | f | F2 hybrid | Both haplotypes recombinant (both parents F1) |

### Model diagram

```
Population A          Population B
     |                     |
     +----------+----------+
                |
               F1 (Mc)
                |
     +----------+----------+
     |                     |
Backcross A (Mb)    Backcross B (Me)
     |                     |
     +----------+----------+
                |
            F2 (Mf)
```

## Recombination Model

For backcross (Mb, Me) and F2 (Mf) models, MONGRAIL accounts for recombination in the F1 parent's gamete. The recombinant haplotype probability is:

```
U(x) = Σ_z U(x|z) · Q(z|d,r)
```

Where:
- `z` is the ancestry vector indicating which loci derive from population A vs B
- `Q(z|d,r)` is the probability of ancestry configuration given marker distances and recombination rate
- `U(x|z)` is the probability of observing haplotype x given ancestry z

The default recombination rate is 1e-8 per base pair (approximately 1 cM per Mb), which can be adjusted with the `-r` option.

## Statistical Framework

MONGRAIL uses a Bayesian approach with:

- **Prior**: Symmetric Dirichlet prior on haplotype frequencies
- **Likelihood**: Posterior predictive distribution (Dirichlet-Multinomial)
- **Inference**: Maximum likelihood classification across models

The posterior predictive probability for observing a haplotype is:

```
P(haplotype | data) = (count + α) / (total + H·α)
```

Where `α = 1/H` is the prior pseudocount and `H` is the number of possible haplotypes.

## Performance

MONGRAIL 2.0 includes several optimizations for fast likelihood computation:

- Pre-computed recombination probabilities Q(z)
- Cached haplotype probabilities U(x)
- Static memory allocation in hot loops
- Compiler optimizations (-O3)

Typical performance: ~0.3 seconds for 31 individuals across 34 chromosomes.

## Examples

### Basic usage

```bash
./mongrail2 examples/simulated_data/popA.GT \
            examples/simulated_data/popB.GT \
            examples/simulated_data/hybrids.GT
```

### With custom recombination rate

```bash
./mongrail2 -r 2e-8 examples/puma/popA.GT \
                    examples/puma/popB.GT \
                    examples/puma/hybrids.GT
```

### Verbose mode with output to file

```bash
./mongrail2 -v popA.GT popB.GT hybrids.GT > results.txt 2> summary.txt
```

## Interpreting Results

1. **Identify the best model**: For each individual, the model with the highest log-likelihood is the most probable classification.

2. **Assess confidence**: Compare log-likelihoods between models. Larger differences indicate stronger evidence.

3. **Consider biological context**:
   - Pure individuals (Ma, Md) should have much higher likelihoods than hybrid models
   - F1 hybrids (Mc) should show intermediate ancestry
   - Backcrosses (Mb, Me) and F2 (Mf) show evidence of recombination

### Example interpretation

```
# indiv      Ma          Mb          Mc          Md          Me          Mf
      0     -0.88       -2.29       -9.02      -12.92      -10.37       -3.64
```

Individual 0: Ma (-0.88) is the best model by a large margin, indicating this individual is most likely pure population A.

## Citation

If you use MONGRAIL in your research, please cite:

> [Citation to be added]

## License

[License to be added]

## Contact

For questions and bug reports, please open an issue on GitHub.

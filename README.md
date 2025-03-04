# EH-STR-parameter-calculator

This script will calculate a number of short tandem repeat (STR) parameters from ExpansionHunter output.

## Input

The vcf and json file resulting from ExpansionHunter version 5. 

## Usage

python EHoutput_parse.py -s $sample -j $sample.json -v $sample.vcf -o $sample.tsv
with python 3.11.8 or higher and python modules json, io, pandas, argparse, logging and math.
Arguments are:
- `-s <arg>` name of the sample
- `-j <arg>` name of the json file that results from ExpansionHunter
- `-v <arg>` name of the vcf file that results from ExpansionHunter
- `-o <arg>` name of the tab separated output file, either .tsv or .txt

## Output

STR parameters extracted from EH vcf and json output are:
- `RepeatID` same as LocusId in EH json file
- `VariantID` same as VariantId in EH json file and REPID in EH vcf file
- `Motif` used in ExpansionHunter, in the reference orientation, same as RepeatUnit in EH json file and RU in EH vcf file
- `ReadLength` average read length in the STR locus, same as ReadLength in the EH json file
- `FragmentLength` average fragment length in the flanking sequence of the STR locus, same as FragmentLength in the EH json file
- `ReadDepth` average read depth in the STR locus, same as Coverage in the EH json file
- `ReadType` this can be 'SPANNING', 'FLANKING' or 'IRR', same as SO in the EH vcf file
- `RepeatNr` repeat number of each allele, same as Genotype in the EH json file and CN in the EH vcf file. The first number is the longest allele, the second is the shortest allele. All other parameters with two numbers follow the corresponding order
- `CI` repeat number confidence interval for each allele, same as GenotypeConfidenceInterval in the EH json file and CI in the EH vcf file
- `SpanningReads` number of spanning reads that is evidence of each allele, same as AD_SP in EH vcf file
- `FlankingReads` number of flanking reads that is evidence of each allele, same as AD_FL in EH vcf file
- `IRRs` number of inrepeat reads that is evidence of each allele, same as AD_IR in EH vcf file

Newly generated STR parameters are:
- `Genotype` this can be 'Hom', 'Het' or 'Hap'
- `AlleleDepth` read depth of each allele calculated from the number of spanning, flanking and inrepeat reads that is evidence for the allele size
- `Missing` whether an allele has not be genotyped (1), zero if not
- `OffDepth` whether an allele has an allele depth that is 5 times higher or lower than the read depth (1), zero if not
- `OffCI` whether an allele has a CI that is larger than the repeat number (1), zero if not
- `FLL` whether the CI of an allele is longer than the fragment length (1), zero if not
- `Flanking` whether an allele is genotyped using flanking reads only (1), zero if not
- `OnlyIRR` whether an allele is genotyped using inrepeat reads (1), zero if not
- `J1C` whether an allele is genotyped using a single consistent read (1), zero if not
- `LCTNC` whether an allele is genotyped with less consistent reads than nonconsistent reads (1), zero if not
- `Qdepth` calculated as the ratio of the allele depth and read depth, number between 0 and 1
- `Qci` calculated as the ratio of the CI and repeat number, number between 0 and 1
- `Qcon` calculated as the ratio of the number of consistent reads and total reads, number between 0 and 1
- `Qnon` calculated as the ratio of the number of nonconsitent reads and total reads, number between 0 and 1

## Publication

This method is described in the following paper:
- Joke J.F.A. van Vugt, Ramona Zwamborn, Egor Dolzhenko and others. The role of disease-associated short tandem repeats in Amyotrophic Lateral Sclerosis, under review.

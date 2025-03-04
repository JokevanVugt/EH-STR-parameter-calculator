# EH-STR-parameter-calculator

This script will calculate a number of short tandem repeat (STR) parameters from ExpansionHunter output.

## Input

The vcf and json file resulting from ExpansionHunter version 5. 

## Usage

python EHoutput_parse.py -s $sample -j $sample.json -v $sample.vcf -o $sample.tsv
with python 3.11.8 or higher and python modules json, io, pandas, argparse, logging and math.
Arguments are:
- -s <arg> name of the sample
- -j <arg> name of the json file that results from ExpansionHunter
- -v <arg> name of the vcf file that results from ExpansionHunter
- -o <arg> name of the tab separated output file, either .tsv or .txt

## Output

The output

## Publication

This method is described in the following paper:
- Joke J.F.A. van Vugt, Ramona Zwamborn, Egor Dolzhenko and others. The role of disease-associated short tandem repeats in Amyotrophic Lateral Sclerosis, under review.

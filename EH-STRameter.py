import json
import io
import pandas as pd
import argparse
import logging
import math

def parse_genotype(genotype, index):
    """Extracts allele information from a genotype string."""
    try:
        return int(genotype.split("/")[index])
    except (IndexError, ValueError):
        return None

def simple_quality(sinfo, motif_length, read_length):
    """Determines the quality of a genotype based on read counts."""
    try:
        flk = eval(sinfo["CountsOfFlankingReads"])
        irr = eval(sinfo["CountsOfInrepeatReads"])
        spn = eval(sinfo["CountsOfSpanningReads"])
        allele1 = parse_genotype(sinfo["Genotype"], 0)
        allele2 = parse_genotype(sinfo["Genotype"], 1)

        if allele1 is not None and allele2 is not None:
            return diploid_quality(flk, irr, spn, allele1, allele2, motif_length, read_length)
        return haploid_quality(flk, irr, spn, allele1, motif_length, read_length)
    
    except (KeyError, SyntaxError):
        return None  # Handle missing or malformed input gracefully

# add consistency numbers to the summaries for comparison 1 and 2 (c1 and c2)
def add_consistency(c1, c2, n1, num, consist_both, consist_just1, consist_just2, nonconsist,consist_other):
    if c1 and c2:
        consist_both += num
    elif c1:
        consist_just1 += num
    elif c2:
        consist_just2 += num
    elif n1:
        nonconsist += num
    else:
        consist_other += num
    return consist_both, consist_just1, consist_just2, nonconsist, consist_other

# add consistency numbers to the summaries for comparison 1 (haploid)
def add_consistency_haploid(c1, n1, num, consist_just1, nonconsist, consist_other):
    if c1:
        consist_just1 += num
    elif n1:
        nonconsist += num
    else:
        consist_other += num
    return consist_just1, nonconsist, consist_other

# identify if flanking reads are consistent with an allele
def flank_consist(flank, allele, motif_length, read_length):
    if flank <= allele:
        return flank * motif_length >= 0.8 * read_length

# identify if flanking reads are not consistent with two alleles
def flank_nonconsist_diploid(flank, allele1, allele2):
    if flank > allele1:
        return flank > allele2

# identify if flanking reads are not consistent with one allele1
def flank_nonconsist(flank, allele):
    return flank > allele

# identify if spanning reads are consistent with an allele
def spn_consist(spn, allele):
    return spn == allele

# identify if spanning reads are not consistent with two alleles
def spn_nonconsist_diploid(spn, allele1, allele2):
    if spn != allele1:
        return spn != allele2

# identify if spanning reads are not consistent with an allele
def spn_nonconsist(spn, allele):
    return spn != allele

# identify if irr's are consistent with an allele
def irr_consist(irr, allele, motif_length, read_length):
    if irr <= allele:
        return irr * motif_length >= 0.8 * read_length

# identify if irr's are not consistent with two alleles
def irr_nonconsist_diploid(irr, allele1, allele2):
    if irr > allele1:
        return irr > allele2

# identify if irr's are not consistent with an allele
def irr_nonconsist(irr, allele):
    return irr > allele

def fix_tuples(test_tuple):
    for test in test_tuple:
        try:
            len_test = len(test)
        except Exception:
            return [test_tuple]
    return test_tuple

# assess the quality of diploid variants based on the number of supporting reads for each allele
def diploid_quality(flk, irr, spn, a1, a2, motif_length, read_length):
    consist_both = consist_just1 = consist_just2 = nonconsist = consist_other = 0

    # fix the tuples
    flk = fix_tuples(flk)
    spn = fix_tuples(spn)
    irr = fix_tuples(irr)

    # categorize the flanking read support
    for f in flk:
        c1 = flank_consist(f[0], a1, motif_length, read_length)
        c2 = flank_consist(f[0], a2, motif_length, read_length)
        n1 = flank_nonconsist_diploid(f[0], a1, a2)
        (consist_both, consist_just1, consist_just2, nonconsist, consist_other,) = add_consistency(c1, c2, n1, f[1], consist_both, consist_just1, consist_just2, nonconsist, consist_other)

    # categorize the spanning read support
    for f in spn:
        c1 = spn_consist(f[0], a1)
        c2 = spn_consist(f[0], a2)
        n1 = spn_nonconsist_diploid(f[0], a1, a2)
        (consist_both, consist_just1, consist_just2, nonconsist, consist_other,) = add_consistency(c1, c2, n1, f[1], consist_both, consist_just1, consist_just2, nonconsist, consist_other)

    # categorize the IRR support
    for f in irr:
        c1 = irr_consist(f[0], a1, motif_length, read_length)
        c2 = irr_consist(f[0], a2, motif_length, read_length)
        n1 = irr_nonconsist_diploid(f[0], a1, a2)
        (consist_both, consist_just1, consist_just2, nonconsist, consist_other,) = add_consistency(c1, c2, n1, f[1], consist_both, consist_just1,consist_just2, nonconsist, consist_other)

    if a1 != a2:
        return_values = [a1, a2, consist_both, consist_just1, consist_just2, nonconsist, consist_other, "Het"]
        return tuple(return_values)

    return_values = [a1, a2, consist_both, consist_just1, consist_just2, nonconsist, consist_other, "Hom"]
    return tuple(return_values)

# assess the quality of haploid variants based on the number of supporting reads for each allele
def haploid_quality(flk, irr, spn, a1, motif_length, read_length):
    consist_just1 = nonconsist = consist_other = 0

    # fix the tuples
    flk = fix_tuples(flk)
    spn = fix_tuples(spn)
    irr = fix_tuples(irr)

    # categorize the flanking read support
    for f in flk:
        c1 = flank_consist(f[0], a1, motif_length, read_length)
        n1 = flank_nonconsist(f[0], a1)
        consist_just1, nonconsist, consist_other = add_consistency_haploid(c1, n1, f[1], consist_just1, nonconsist, consist_other)

    # categorize the spanning read support
    for f in spn:
        c1 = spn_consist(f[0], a1)
        n1 = spn_nonconsist(f[0], a1)
        consist_just1, nonconsist, consist_other = add_consistency_haploid(c1, n1, f[1], consist_just1, nonconsist, consist_other)

    # categorize the IRR support
    for f in irr:
        c1 = irr_consist(f[0], a1, motif_length, read_length)
        n1 = irr_nonconsist(f[0], a1)
        consist_just1, nonconsist, consist_other = add_consistency_haploid(c1, n1, f[1], consist_just1, nonconsist, consist_other)

    return a1, consist_just1, nonconsist, consist_other, "Hap"

def read_vcf(sample, vcffile, eh_sampleid):
    with open(vcffile) as f:
        lines = [l for l in f if not l.startswith("##")]

    eh_input = pd.read_csv(
        io.StringIO("".join(lines)), sep="\t",
        dtype={"#CHROM": str, "POS": int, "ID": str, "REF": str,
               "ALT": str, "QUAL": str, "FILTER": str, "INFO": str}
    ).rename(columns={"#CHROM": "CHROM"})

    # Extract INFO fields
    info_fields = ["END", "REF", "RL", "RU", "VariantID", "RepeatID"]
    eh_input[info_fields] = eh_input["INFO"].astype(str).str.split(";", expand=True)

    # Extract sample-specific fields
    sample_fields = ["GT", "ReadType", "RepeatNr", "CI", "SpanningReads", "FlankingReads", "IRRs", "LC"]
    eh_input[sample_fields] = eh_input[eh_sampleid].astype(str).str.split(":", expand=True)

    # Clean up VariantID and RepeatID
    eh_input["VariantID"] = eh_input["VariantID"].apply(lambda x: x.replace("VARID=", "") if isinstance(x, str) else x)
    eh_input["RepeatID"] = eh_input["RepeatID"].apply(lambda x: x.replace("REPID=", "") if isinstance(x, str) else x)

    eh_input["SampleID"] = sample

    return eh_input[["VariantID", "RepeatNr", "CI", "ReadType", "SpanningReads", "FlankingReads", "IRRs"]]

def split_values(line_info):
    """Splits the input strings into lists of values"""
    return {
        "RepeatNr": [int(x) for x in line_info["RepeatNr"].split("/")],
        "CI": [x.split("-") for x in line_info["CI"].split("/")],
        "ReadType": line_info["ReadType"].split("/"),
        "SpanningReads": [int(x) for x in line_info["SpanningReads"].split("/")],
        "FlankingReads": [int(x) for x in line_info["FlankingReads"].split("/")],
        "IRRs": [int(x) for x in line_info["IRRs"].split("/")]
    }

def calculate_allele_depth(RepeatNr, ReadType, SpanningReads, FlankingReads, IRRs, Genotype, Motif):
    """Calculates allele depth based on ReadType and Genotype"""
    RepeatNr_divide = RepeatNr if RepeatNr != 0 else 0.01
    if ReadType == "SPANNING":
        depth = (SpanningReads + (FlankingReads / 4)) / (2 if Genotype == "Hom" else 1)
    elif ReadType == "INREPEAT":
        depth = ((IRRs + (FlankingReads / 4)) * (75 if Genotype == "Hom" else 150)) / (RepeatNr_divide * len(Motif))
    elif ReadType == "FLANKING":
        depth = FlankingReads / (8 if Genotype == "Hom" else 4)
    else:
        return "NA"
    return round(depth, 2)

def adjust_for_ATXN1(value, variant_id):
    """Applies ATXN1 correction if necessary."""
    if value == "NA":
        return "NA"
    return value - 1 if variant_id == "ATXN1" else value

def setGenotype_to_NA(swap):
    """Handles NA Genotype by setting all values to NA"""
    swap.update({
        "RepeatNr_max": "NA", "RepeatNr_min": "NA", "CI_max_1": "NA", "CI_max_2": "NA",
        "CI_min_1": "NA", "CI_min_2": "NA", "ReadType_max": "NA", "ReadType_min": "NA",
        "SpanningReads_max": "NA", "SpanningReads_min": "NA", "FlankingReads_max": "NA",
        "FlankingReads_min": "NA", "IRRs_max": "NA", "IRRs_min": "NA", "OnlyIRR": "NA",
        "Consist_just_max": "NA", "Consist_just_min": "NA", "Consist_both": "NA", 
        "Consist_other": "NA", "Nonconsist": "NA", "AlleleDepth_max": "NA",
        "AlleleDepth_min": "NA", "Missing": 2, "Missing_max": 1, "Missing_min": 1,
        "OffAllelicDepth_max": "NA", "OffAllelicDepth_min": "NA", "OffAllelicDepth": "NA",
        "FLL_max": "NA", "FLL_min": "NA", "FragmentLengthLimited": "NA",
        "Flanking_max": "NA", "Flanking_min": "NA", "IfFlanking": "NA",
        "J1C_max": "NA", "J1C_min": "NA", "J1C": "NA",
        "LCTNC_max": "NA", "LCTNC_min": "NA", "LCTNC": "NA",
        "CIratio_max": "NA", "CIratio_min": "NA", "CIratio": "NA",
        "Qcon_max": "NA", "Qcon_min": "NA", "Qnon": "NA", "Qci_max": "NA",
        "Qci_min": "NA", "Qdepth_max": "NA", "Qdepth_min": "NA"
    })

def calculate_missing(RepeatNr_max, RepeatNr_min, Genotype):
    """Returns the sum of missing values (1 if missing, 0 otherwise), 
    considering RepeatNr_min only if Genotype is 'Het' or 'Hom'
    """
    missing_max = int(RepeatNr_max == "NA" or RepeatNr_max == ".")
    missing_min = int((RepeatNr_min == "NA" or RepeatNr_min == ".") and Genotype in ["Het", "Hom"])
    return missing_max + missing_min

def calculate_off_allelic_depth(allele_depth, read_depth, missing):
    """Returns 1 if OffAllelicDepth condition is met, otherwise 0"""
    if allele_depth in [None, "NA"]:
        allele_depth = 0
    return int(
        ((read_depth + 1) / (allele_depth + 1) > 5 and missing == 0) or
        ((read_depth + 1) / (allele_depth + 1) < 0.2 and missing == 0) or
        (allele_depth == 0 and missing == 0)
    )

def calculate_fragment_length_limited(fragment_length, CI2, missing, motif):
    """Returns 1 if CI_2 >= fragment length, 0 otherwise"""
    if CI2 in [None, "NA"]:
        CI2 == 0
    return int(int(CI2)*len(motif) >= float(fragment_length) and missing == 0)

def calculate_if_flanking(readtype, missing):
    """Returns 1 if ReadType == FLANKING, 0 otherwise"""
    if readtype in [None, "NA"]:
        readtype == 0
    return int(readtype == "FLANKING" and missing == 0)

def calculate_J1C(consist_both, consist_one, missing, genotype):
    """Returns 1 if just one consistent read, otherwise 0"""
    if consist_both in [None, "NA"]:
        consist_both == 0
    return int(
        (genotype == "Hom" and consist_both/2 < 2 and missing == 0) or
        (genotype == "Het" and consist_one + consist_both/2 < 2 and missing == 0) or
        (genotype == "Hap" and consist_one < 2 & missing == 0)
    )

def calculate_LCTNC(consist_both, consist_one, nonconsist, missing, genotype):
    """Returns 1 if less consistent than nonconsistent reads, otherwise 0"""
    if nonconsist in [None, "NA"]:
        nonconsist == 0
    return int(
        (genotype == "Hom" and consist_both/2 < nonconsist and missing == 0) or
        (genotype == "Het" and consist_one + consist_both/2 < nonconsist and missing == 0) or
        (genotype == "Hap" and consist_one < nonconsist & missing == 0)
    )

def calculate_CIratio(repeatnr, CI_1, CI_2, missing):
    """Returns 1 if less consistent than nonconsistent reads, otherwise 0"""
    if repeatnr in [None, "NA"]:
        CI_2 == 1 and CI_1 == 1 and repeatnr == 1
    return int((int(CI_2) - int(CI_1)) > int(repeatnr) and missing == 0)

def merge_json_vcf(line_info):
    """Merges json and vcf info and returns structured data."""
    swap = {}

    if line_info["Genotype"] == "NA":
        setGenotype_to_NA(swap)
        return swap

    values = split_values(line_info)
    RepeatNr = values["RepeatNr"]
    
    if len(RepeatNr) == 1:
        swap.update({
            "Genotype": line_info["Genotype"], "FragmentLength": line_info["FragmentLength"], "Motif": line_info["Motif"],
            "ReadDepth": line_info["ReadDepth"], "RepeatNr_max": RepeatNr[0], "RepeatNr_min": "NA",
            "CI_max_1": values["CI"][0][0], "CI_max_2": values["CI"][0][1],
            "CI_min_1": "NA", "CI_min_2": "NA",
            "ReadType_max": values["ReadType"][0], "ReadType_min": "NA",
            "SpanningReads_max": values["SpanningReads"][0], "SpanningReads_min": "NA",
            "FlankingReads_max": values["FlankingReads"][0], "FlankingReads_min": "NA",
            "IRRs_max": values["IRRs"][0], "IRRs_min": "NA",
            "OnlyIRR": int(values["ReadType"][0] == "INREPEAT"),
            "Consist_just_max": "NA" if line_info["Allele1"] == line_info["Allele2"] == "NA" else line_info["Consist_just1"],
            "Consist_just_min": "NA", "Consist_both": line_info["Consist_both"], "Nonconsist": line_info["Nonconsist"],
            "Consist_other": line_info["Consist_other"]
        })
        
        swap["AlleleDepth_max"] = calculate_allele_depth(
            swap["RepeatNr_max"], swap["ReadType_max"],
            swap["SpanningReads_max"], swap["FlankingReads_max"],
            swap["IRRs_max"], line_info["Genotype"], line_info["Motif"]
        )
        swap["AlleleDepth_min"] = "NA"

    elif len(RepeatNr) == 2:
        RepeatNr_max = max(RepeatNr)
        max_idx = RepeatNr.index(RepeatNr_max)
        min_idx = 1 - max_idx  # Opposite index

        swap.update({
            "Genotype": line_info["Genotype"], "FragmentLength": line_info["FragmentLength"], "Motif": line_info["Motif"],
            "ReadDepth": line_info["ReadDepth"], "RepeatNr_max": RepeatNr_max, "RepeatNr_min": RepeatNr[min_idx],
            "CI_max_1": int(values["CI"][max_idx][0]), "CI_max_2": int(values["CI"][max_idx][1]),
            "CI_min_1": int(values["CI"][min_idx][0]), "CI_min_2": int(values["CI"][min_idx][1]),
            "ReadType_max": values["ReadType"][max_idx], "ReadType_min": values["ReadType"][min_idx],
            "SpanningReads_max": values["SpanningReads"][max_idx], "SpanningReads_min": values["SpanningReads"][min_idx],
            "FlankingReads_max": values["FlankingReads"][max_idx], "FlankingReads_min": values["FlankingReads"][min_idx],
            "IRRs_max": values["IRRs"][max_idx], "IRRs_min": values["IRRs"][min_idx],
            "OnlyIRR": int(values["ReadType"][max_idx] == "INREPEAT" and values["ReadType"][min_idx] == "INREPEAT"),
            "Consist_both": line_info["Consist_both"], "Nonconsist": line_info["Nonconsist"], "Consist_other": line_info["Consist_other"]
        })

        if line_info["Allele1"] == "NA" == line_info["Allele2"]:
            swap["Consist_just_max"] = "NA"
            swap["Consist_just_min"] = "NA"
        elif line_info["Allele2"] == "NA":
            swap["Consist_just_max"] = line_info["Consist_just1"]
            swap["Consist_just_min"] = "NA"
        elif line_info["Allele2"] == RepeatNr_max:
            swap["Consist_just_max"] = line_info["Consist_just2"]
            swap["Consist_just_min"] = line_info["Consist_just1"]
        else:
            swap["Consist_just_max"] = line_info["Consist_just1"]
            swap["Consist_just_min"] = line_info["Consist_just2"]

        swap["AlleleDepth_max"] = calculate_allele_depth(
            swap["RepeatNr_max"], swap["ReadType_max"],
            swap["SpanningReads_max"], swap["FlankingReads_max"],
            swap["IRRs_max"], line_info["Genotype"], line_info["Motif"]
        )

        swap["AlleleDepth_min"] = calculate_allele_depth(
            swap["RepeatNr_min"], swap["ReadType_min"],
            swap["SpanningReads_min"], swap["FlankingReads_min"],
            swap["IRRs_min"], line_info["Genotype"], line_info["Motif"]
        )

    else:
        setGenotype_to_NA(swap)
        return swap

    # ATXN1 correction
    swap["RepeatNr_max"] = adjust_for_ATXN1(swap["RepeatNr_max"], line_info["VariantID"])
    swap["CI_max_1"] = adjust_for_ATXN1(swap["CI_max_1"], line_info["VariantID"])
    swap["CI_max_2"] = adjust_for_ATXN1(swap["CI_max_2"], line_info["VariantID"])
    swap["RepeatNr_min"] = adjust_for_ATXN1(swap["RepeatNr_min"], line_info["VariantID"])
    swap["CI_min_1"] = adjust_for_ATXN1(swap["CI_min_1"], line_info["VariantID"])
    swap["CI_min_2"] = adjust_for_ATXN1(swap["CI_min_2"], line_info["VariantID"])

    # Missing
    swap["Missing"] = calculate_missing(swap["RepeatNr_max"], swap["RepeatNr_min"], swap["Genotype"])
    swap["Missing_max"] = int(swap["RepeatNr_max"] == "NA")
    swap["Missing_min"] = int(swap["RepeatNr_min"] == "NA" and swap["Genotype"] in ["Het", "Hom"])

    # OffAllelicDepth
    ReadDepth = line_info.get("ReadDepth", 0)
    swap["OffAllelicDepth_max"] = calculate_off_allelic_depth(swap["AlleleDepth_max"], ReadDepth, swap["Missing_max"])
    swap["OffAllelicDepth_min"] = 0 if swap["Genotype"] == "Hap" else calculate_off_allelic_depth(swap["AlleleDepth_min"], ReadDepth, swap["Missing_min"])
    swap["OffAllelicDepth"] = min(2, swap["OffAllelicDepth_max"] + swap["OffAllelicDepth_min"])

    # FragmentLengthLimited
    swap["FLL_max"] = calculate_fragment_length_limited(swap["FragmentLength"], swap["CI_max_2"], swap["Missing_max"], swap["Motif"])
    swap["FLL_min"] = 0 if swap["Genotype"] == "Hap" else calculate_fragment_length_limited(swap["FragmentLength"], swap["CI_min_2"], swap["Missing_min"], swap["Motif"])
    swap["FragmentLengthLimited"] = min(2, swap["FLL_max"] + swap["FLL_min"])

    # IfFlanking
    swap["Flanking_max"] = calculate_if_flanking(swap["ReadType_max"], swap["Missing_max"])
    swap["Flanking_min"] = 0 if swap["Genotype"] == "Hap" else calculate_if_flanking(swap["ReadType_min"], swap["Missing_min"])
    swap["IfFlanking"] = min(2, swap["Flanking_max"] + swap["Flanking_min"])

    # J1C
    swap["J1C_max"] = calculate_J1C(swap["Consist_both"], swap["Consist_just_max"], swap["Missing_max"], swap["Genotype"])
    swap["J1C_min"] = 0 if swap["Genotype"] == "Hap" else calculate_J1C(swap["Consist_both"], swap["Consist_just_min"], swap["Missing_min"], swap["Genotype"])
    swap["J1C"] = min(2, swap["J1C_max"] + swap["J1C_min"])

    # LCTNC
    swap["LCTNC_max"] = calculate_LCTNC(swap["Consist_both"], swap["Consist_just_max"], swap["Nonconsist"], swap["Missing_max"], swap["Genotype"])
    swap["LCTNC_min"] = 0 if swap["Genotype"] == "Hap" else calculate_LCTNC(swap["Consist_both"], swap["Consist_just_min"], swap["Nonconsist"], swap["Missing_min"], swap["Genotype"])
    swap["LCTNC"] = min(2, swap["LCTNC_max"] + swap["LCTNC_min"])

    # deltaCI larger than RepeatNr
    swap["CIratio_max"] = calculate_CIratio(swap["RepeatNr_max"], swap["CI_max_1"], swap["CI_max_2"], swap["Missing_max"])
    swap["CIratio_min"] = 0 if swap["Genotype"] == "Hap" else calculate_CIratio(swap["RepeatNr_min"], swap["CI_min_1"], swap["CI_min_2"], swap["Missing_min"])
    swap["CIratio"] = min(2, swap["CIratio_max"] + swap["CIratio_min"])

    # Qcon
    swap["Consist_both_half"] = swap["Consist_both"]/2 if swap["Consist_both"] != 'NA' else 0
    if swap["Genotype"] == "Hom":
        swap["Qcon_max"] = round(swap["Consist_both"] / (swap["SpanningReads_max"] + swap["FlankingReads_max"] + swap["IRRs_max"]), 2)
        swap["Qcon_min"] = round(swap["Consist_both"] / (swap["SpanningReads_min"] + swap["FlankingReads_min"] + swap["IRRs_min"]), 2)
    elif swap["Genotype"] == "Het":
        swap["Qcon_max"] = round((swap["Consist_both_half"] + swap["Consist_just_max"])/(swap["SpanningReads_max"] + swap["FlankingReads_max"] + swap["IRRs_max"]), 2)
        swap["Qcon_min"] = round((swap["Consist_both_half"] + swap["Consist_just_min"])/(swap["SpanningReads_min"] + swap["FlankingReads_min"] + swap["IRRs_min"]), 2)
    elif swap["Genotype"] == "Hap":
        swap["Qcon_max"] = round(swap["Consist_just_max"]/(swap["SpanningReads_max"] + swap["FlankingReads_max"] + swap["IRRs_max"]), 2)
        swap["Qcon_min"] = 'NA'
    else:
        swap["Qcon_max"] = 'NA'
        swap["Qcon_min"] = 'NA'

    # Qnon
    if swap["Nonconsist"] != 'NA':
        list = [swap["Consist_both"], swap["Consist_just_max"], swap["Consist_just_min"], swap["Nonconsist"], swap["Consist_other"]]
        filtered_list = [x for x in list if x != 'NA']
        swap["Qnon"] = round(1/math.exp(4*swap["Nonconsist"]/sum(filtered_list)), 2)
    else:
        swap["Qnon"] = 'NA'

    # Qci
    if swap["RepeatNr_max"] != 'NA':
        swap["Qci_max"] = round(1/math.exp(4*(int(swap["CI_max_2"]) - int(swap["CI_max_1"]))/swap["RepeatNr_max"]), 2)
    else:
        swap["Qci_max"] = 'NA'
    if swap["RepeatNr_min"] != 'NA':
        swap["Qci_min"] = round(1/math.exp(4*(int(swap["CI_min_2"]) - int(swap["CI_min_1"]))/swap["RepeatNr_min"]), 2)
    else:
        swap["Qci_min"] = 'NA'

    # Qdepth
    if swap["AlleleDepth_max"] != 'NA':
        if swap["AlleleDepth_max"] >= swap["ReadDepth"]:
            swap["Qdepth_max"] = round(swap["ReadDepth"]/swap["AlleleDepth_max"] ,2)
        else:
            swap["Qdepth_max"] = round(swap["AlleleDepth_max"]/swap["ReadDepth"] ,2)
    else:
        swap["Qdepth_max"] = 'NA'
    if swap["AlleleDepth_min"] != 'NA':
        if swap["AlleleDepth_min"] >= swap["ReadDepth"]:
            swap["Qdepth_min"] = round(swap["ReadDepth"]/swap["AlleleDepth_min"] ,2)
        else:
            swap["Qdepth_min"] = round(swap["AlleleDepth_min"]/swap["ReadDepth"] ,2)    
    else:
        swap["Qdepth_min"] = 'NA'

    return swap

logging.basicConfig(format="%(asctime)s %(message)s", level=logging.DEBUG, datefmt="%Y-%m-%d %H:%M:%S")

def main():
    parser = argparse.ArgumentParser(
        description="Parses EH json and vcf into single txt file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-s", "--samplename", help="This is the sample name")
    parser.add_argument("-j", "--jsonfile", help="This is the json file from EH")
    parser.add_argument("-v", "--vcffile", help="This is the vcf file from EH")
    parser.add_argument("-o", "--outputfile", help="This is the output file")
    args = parser.parse_args()

    sample, jsonfile, vcffile, outputfile = args.samplename, args.jsonfile, args.vcffile, args.outputfile

#    logging.info("Loading json file")
    with open(jsonfile) as json_input_file:
        json_input = json.load(json_input_file)
    eh_sampleid = json_input["SampleParameters"]["SampleId"]
    json_input = json_input["LocusResults"]

#    logging.info("Loading VCF file")
    vcf_df = read_vcf(sample, vcffile, eh_sampleid)
    vcf_df.set_index(["VariantID"], inplace=True)

#    logging.info("Looping over JSON items")
    i = 0
    fh_out = open(outputfile, "w")
    # write header
    fh_out.write(
        "SampleID\tRepeatID\tVariantID\tMotif\tReadLength"
        "\tFragmentLength\tReadDepth\tGenotype\tReadType"
        "\tRepeatNr\tCI\tSpanningReads\tFlankingReads\tIRRs"
        "\tAlleleDepth\tMissing\tOffDepth\tOffCI\tFLL"
        "\tFlanking\tOnlyIRR\tJ1C\tLCTNC\tQdepth\tQci\tQcon\tQnon\n"
    )

    no_genotype_dict = {key: "NA" for key in [
        "Genotype", "Allele1", "Allele2", "Consist_both", "Consist_just1", "Consist_just2", "Nonconsist", "Consist_other"
    ]}

    for _, locus_rec in json_input.items():
        i = i + 1

        if "Variants" not in locus_rec:
            logging.warning(f'No variants in {locus_rec["LocusId"]} of {sample}')
            continue

        line_info = {
            "SampleID": sample,
            "ReadDepth": locus_rec["Coverage"],
            "FragmentLength": locus_rec["FragmentLength"],
            "ReadLength": locus_rec["ReadLength"],
            "RepeatID": locus_rec["LocusId"],
        }

        if locus_rec["AlleleCount"] == 2 and locus_rec["ReadLength"] > 0:
            for variant_id, variant_rec in locus_rec["Variants"].items():
                line_info["VariantID"] = variant_id
                motif = variant_rec["RepeatUnit"]
                motif_length = len(motif)
                line_info["Motif"] = motif

                if "Genotype" in variant_rec:
                    read_consist = simple_quality(variant_rec, motif_length, locus_rec["ReadLength"])
                    line_info["Allele1"] = read_consist[0]
                    line_info["Allele2"] = read_consist[1]
                    line_info["Consist_both"] = read_consist[2]
                    line_info["Consist_just1"] = read_consist[3]
                    line_info["Consist_just2"] = read_consist[4]
                    line_info["Nonconsist"] = read_consist[5]
                    line_info["Consist_other"] = read_consist[6]
                    line_info["Genotype"] = read_consist[7]

                else:
                    line_info.update(no_genotype_dict)

                line_info.update(vcf_df.loc[variant_id].to_dict())
                line_info.update(merge_json_vcf(line_info))

                # write to output
                try:
                    fh_out.write(
                        f'{line_info["SampleID"]}\t{line_info["RepeatID"]}\t{line_info["VariantID"]}\t{line_info["Motif"]}\t{line_info["ReadLength"]}'
                        f'\t{line_info["FragmentLength"]}\t{line_info["ReadDepth"]:.2f}\t{line_info["Genotype"]}\t{line_info["ReadType_max"]}/{line_info["ReadType_min"]}'
                        f'\t{line_info["RepeatNr_max"]}/{line_info["RepeatNr_min"]}\t{line_info["CI_max_1"]}-{line_info["CI_max_2"]}/{line_info["CI_min_1"]}-{line_info["CI_min_2"]}'
                        f'\t{line_info["SpanningReads_max"]}/{line_info["SpanningReads_min"]}\t{line_info["FlankingReads_max"]}/{line_info["FlankingReads_min"]}'
                        f'\t{line_info["IRRs_max"]}/{line_info["IRRs_min"]}\t{line_info["AlleleDepth_max"]}/{line_info["AlleleDepth_min"]}\t{line_info["Missing"]}'
                        f'\t{line_info["OffAllelicDepth"]}\t{line_info["CIratio"]}\t{line_info["FragmentLengthLimited"]}\t{line_info["IfFlanking"]}\t{line_info["OnlyIRR"]}'
                        f'\t{line_info["J1C"]}\t{line_info["LCTNC"]}\t{line_info["Qdepth_max"]}/{line_info["Qdepth_min"]}\t{line_info["Qci_max"]}/{line_info["Qci_min"]}'
                        f'\t{line_info["Qcon_max"]}/{line_info["Qcon_min"]}\t{line_info["Qnon"]}\n'
                    )
                except KeyError as e:
                    logging.warning(e)
                    logging.warning(line_info)
                    raise

        elif locus_rec["AlleleCount"] == 1 and locus_rec["ReadLength"] > 0:
            for variant_id, variant_rec in locus_rec["Variants"].items():
                line_info["VariantID"] = variant_id
                motif = variant_rec["RepeatUnit"]
                motif_length = len(motif)
                line_info["Motif"] = motif

                if "Genotype" in variant_rec:
                    repeat_len = int(variant_rec["Genotype"])
                    CI_max = variant_rec["GenotypeConfidenceInterval"]
                    CI_max_1 = int(CI_max.split("-")[0])
                    CI_max_2 = int(CI_max.split("-")[1])
                    line_info["RepeatNr_max"] = repeat_len
                    line_info["CI_max_1"] = CI_max_1
                    line_info["CI_max_2"] = CI_max_2
                    read_consist = simple_quality(variant_rec, motif_length, locus_rec["ReadLength"])
                    line_info["Allele1"] = read_consist[0]
                    line_info["Consist_just1"] = read_consist[1]
                    line_info["Nonconsist"] = read_consist[2]
                    line_info["Consist_other"] = read_consist[3]
                    line_info["Genotype"] = read_consist[4]
                    # following fields are set to NA since there is only one allel
                    line_info["Consist_both"] = "NA"
                    line_info["Consist_just2"] = "NA"
                    line_info["Allele2"] = "NA"
                    line_info["CI_min_1"] = "NA"
                    line_info["CI_min_2"] = "NA"

                else:
                    line_info.update(no_genotype_dict)

                line_info.update(vcf_df.loc[variant_id].to_dict())
                line_info.update(merge_json_vcf(line_info))

                # write to output
                try:
                    fh_out.write(
                        f'{line_info["SampleID"]}\t{line_info["RepeatID"]}\t{line_info["VariantID"]}\t{line_info["Motif"]}\t{line_info["ReadLength"]}'
                        f'\t{line_info["FragmentLength"]}\t{line_info["ReadDepth"]:.2f}\t{line_info["Genotype"]}\t{line_info["ReadType_max"]}'
                        f'\t{line_info["RepeatNr_max"]}\t{line_info["CI_max_1"]}-{line_info["CI_max_2"]}'
                        f'\t{line_info["SpanningReads_max"]}\t{line_info["FlankingReads_max"]}\t{line_info["IRRs_max"]}'
                        f'\t{line_info["AlleleDepth_max"]}\t{line_info["Missing"]}\t{line_info["OffAllelicDepth"]}\t{line_info["CIratio"]}'
                        f'\t{line_info["FragmentLengthLimited"]}\t{line_info["IfFlanking"]}\t{line_info["OnlyIRR"]}\t{line_info["J1C"]}\t{line_info["LCTNC"]}'
                        f'\t{line_info["Qdepth_max"]}\t{line_info["Qci_max"]}'
                        f'\t{line_info["Qcon_max"]}\t{line_info["Qnon"]}\n'
                    )
                except KeyError as e:
                    logging.warning(e)
                    logging.warning(line_info)
                    raise

        else:
            logging.warning(f'No genotyping of {locus_rec["LocusId"]} of {sample}')
            continue

    fh_out.close()
    logging.info(f"{sample} done")

if __name__ == "__main__":
    main()

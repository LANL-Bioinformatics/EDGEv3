#!/usr/bin/env python

import json
import argparse
import re


def get_args():
    parser = argparse.ArgumentParser(description = "Convert FaQCs stats report to JSON")
    parser.add_argument("-i", "--input", help="Input file", required = True)
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()
    in_file = args.input
    json_dict = {}
    with open(in_file, 'rt') as f:
        after_trimming = False #disambiguate lines with the same starting text
        for line in f:
            line = line.strip()
            if line == "After Trimming":
                after_trimming = True
            m = re.search(r': ([\d.]+)', line)
            if m:
                if line.startswith("Reads #") and not after_trimming:
                    json_dict["inputReads"] = m.group(1)
                elif line.startswith("Total bases") and not after_trimming:
                    json_dict["inputBases"] = m.group(1)
                elif line.startswith("Reads Length") and not after_trimming:
                    json_dict["inputReadLength"] = m.group(1)
                elif line.startswith("Reads #") and after_trimming:
                    json_dict["outputReads"] = m.group(1)
                elif line.startswith("Total bases") and after_trimming:
                    json_dict["outputBases"] = m.group(1)
                elif line.startswith("Mean Reads Length") and after_trimming:
                    json_dict["outputReadLength"] = m.group(1)
                elif line.startswith("Paired Reads #"):
                    json_dict["outputPairedReads"] = m.group(1)
                elif line.startswith("Paired total bases"):
                    json_dict["outputPairedBases"] = m.group(1)
                elif line.startswith("Unpaired Reads #"):
                    json_dict["outputUnpairedReads"] = m.group(1)
                elif line.startswith("Unpaired total bases"):
                    json_dict["outputUnpairedBases"] = m.group(1)
                elif line.startswith("Discarded reads #"):
                    json_dict["filteredTotalReads"] = m.group(1)
                elif line.startswith("Trimmed bases"):
                    json_dict["trimmedTotalBases"] = m.group(1)
                elif line.startswith("Reads Filtered by length"):
                    json_dict["lenFilteredReads"] = m.group(1)
                elif line.startswith("Bases Filtered by length"):
                    json_dict["lenFilteredBases"] = m.group(1)
                elif line.startswith("Reads Filtered by continuous"):
                    json_dict["nFilteredReads"] = m.group(1)
                elif line.startswith("Bases Filtered by continuous"):
                    json_dict["nFilteredBases"] = m.group(1)
                elif line.startswith("Reads Filtered by low complexity"):
                    json_dict["lcFilteredReads"] = m.group(1)
                elif line.startswith("Bases Filtered by low complexity"):
                    json_dict["lcFilteredBases"] = m.group(1)
                elif line.startswith("Reads Trimmed by quality"):
                    json_dict["qualTrimmedReads"] = m.group(1)
                elif line.startswith("Bases Trimmed by quality"):
                    json_dict["qualTrimmedBases"] = m.group(1)
    with open("QC.stats.json", 'w') as f:
        json.dump(json_dict, f)


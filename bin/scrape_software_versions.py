#!/usr/bin/env python

from __future__ import print_function
from collections import OrderedDict
import re

# add regexes for tools in process get_software_versions
regexes = {
    "nf-core/rvfvtyping": ["pipeline.version.txt", r"(\S+)"],
    "Nextflow": ["nextflow.version.txt", r"(\S+)"],
    "ncbi-genome-download":  ["ncbi-genome-download.version.txt", r"(\S+)"],
    "DIAMOND":  ["diamond.version.txt",r"(\S+)"],
    "MAFFT":    ["mafft.version.txt", r"v(\S+)"],
    "IQTREE":   ["iqtree.version.txt", r"(\S+)"],
    "biopython": ["biopython.version.txt", r"(\S+)"],
    "dendropy":  ["dendropy.version.txt", r"(\S+)"],
    "seqmagick":  ["seqmagick.version.txt",r"v(\S+)"],
}
    #"MultiQC": ["v_multiqc.txt", r"multiqc, version (\S+)"],
results = OrderedDict()
results["nf-core/rvfvtyping"] = '<span style="color:#999999;">N/A</span>'
results["Nextflow"] = '<span style="color:#999999;">N/A</span>'
results["ncbi-genome-download"] = '<span style="color:#999999;">N/A</span>'
results["DIAMOND"] = '<span style="color:#999999;">N/A</span>'
results["MAFFT"] = '<span style="color:#999999;">N/A</span>'
results["IQTREE"] = '<span style="color:#999999;">N/A</span>'
results["biopython"] = '<span style="color:#999999;">N/A</span>'
results["dendropy"] = '<span style="color:#999999;">N/A</span>'
results["seqmagick"] = '<span style="color:#999999;">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/rvfvtyping Software Versions'
section_href: 'https://github.com/nf-core/rvfvtyping'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))

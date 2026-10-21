#!/usr/bin/env python3
"""Write <outdir>/versions.tsv: env, package, version, build, channel for every
environment under /opt/envs (conda packages, plus the pip-installed crispro and
protiler as micromamba reports them, plus the R package breakfast).
Run at image build time (see Dockerfile); not needed at run time."""
import json
import os
import subprocess
import sys

ENVS = ["r", "py", "ngs", "mageck", "crispro", "protiler"]
out_dir = sys.argv[1]

rows = []
for env in ENVS:
    prefix = "/opt/envs/" + env
    raw = subprocess.run(["micromamba", "list", "-p", prefix, "--json"],
                         check=True, capture_output=True, text=True).stdout
    for p in sorted(json.loads(raw), key=lambda d: d["name"]):
        channel = p.get("channel") or p.get("base_url") or ""
        rows.append((env, p["name"], p["version"], p.get("build_string", ""), channel))

# R packages that are not conda packages
rows.append(("protiler", "breakfast", "0.1.0", "R CMD INSTALL", "CRAN archive"))

with open(os.path.join(out_dir, "versions.tsv"), "w") as fh:
    fh.write("env\tpackage\tversion\tbuild\tchannel\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
print("versions.tsv: %d rows" % len(rows))

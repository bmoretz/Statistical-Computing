#First add depencies for the example
using Pkg; Pkg.add.(["Plots", "DSP"])
using Weave

list_out_formats()

weave("FIR_design.md", doctype=:"github", out_path=:pwd)

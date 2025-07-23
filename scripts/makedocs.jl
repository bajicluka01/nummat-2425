using Weave

#Weave.weave("scripts/demo.jl", doctype="minted2pdf", out_path="report/report.pdf")
Weave.weave("scripts/demo.jl", doctype="md2pdf", out_path="report/report.pdf")
#Weave.weave("report/report.tex", doctype="minted2pdf", out_path="report/report.pdf")

using Weave;

list_out_formats()

weave("01_Introduciton.jmd", doctype=:"pandoc2pdf", out_path=:pwd)
weave("02_SimpleLinearOptimization.jmd", doctype=:"pandoc2pdf", out_path=:pwd)
weave("03_Julia-Basics.jmd", doctype=:"pandoc2pdf", out_path=:pwd)

install.packages('tinytex')

tinytex::install_tinytex()

# tinytex:::is_tinytex()

library(tinytex)
tlmgr_search('/times.sty')   # search for times.sty
tlmgr_install('psnfss')      # install the psnfss package
tlmgr_update()               # update everything


tinytex::install_tinytex()

tinytex::uninstall_tinytex()

# LIBRARY-------
library("meta", "metasens")
settings.meta(digits=2)
# IMPORT DATA
joy = read.csv(here::here('test_meta.txt')) 
str(joy)
# resp h = medoc
# resp p = placebo
# fail = failure to placebo if p, to medoc if p
joy$miss = ifelse((joy$drop.h + joy$drop.p) == 0, 
                    c("Without missing data"), c("With missing data"))

m.publ = metabin (resp.h, resp.h + fail.h, resp.p, resp.p + fail.p, 
                  data = joy, studlab = paste0(author,"(", year, ")"), 
                  method.tau = "PM") 
m.publ

forest(m.publ, sortvar=year, prediction=TRUE, label.left =
          "Favours placebo", label.right = "Favours haloperidol")

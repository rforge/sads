(TeX-add-style-hook "sads_tutorial"
 (lambda ()
    (LaTeX-add-labels
     "eq:poisson"
     "eq:integral"
     "equation"
     "eq:poiexp"
     "fig:fig1")
    (TeX-add-symbols
     '("code" 1)
     "R")
    (TeX-run-style-hooks
     "graphicx"
     "fontenc"
     "T1"
     "inputenc"
     "utf8"
     "babel"
     "brazil"
     "latex2e"
     "art10"
     "article")))


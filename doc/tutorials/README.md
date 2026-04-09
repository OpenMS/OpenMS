# Tutorials
Handouts and workflows as used in the tutorial session during user meetings. 

Build handout:
- compile PDF from `Handout/handout.tex` using XeLatex and BibTex (for citations) with `Handout` as working directory:
  `cd Handout && xelatex -shell-escape handout && bibtex handout && xelatex -shell-escape handout && xelatex -shell-escape handout`


Note:
- Please put only the handouts and workflows (small files) here. Installers, external tools, and example data should be kept in the releases.
- Please assign a Git tag for any special tutorial event with a separate release

.PHONY: all

all:
	rm -f *.bbl
	cd manuscript && latexmk -gg -pdf main.tex
	cd manuscript && latexmk -c

%.png: %.pdf
	convert -flatten $< $@

manuscript/figures/prospect_pairs_N.png: manuscript/figures/prospect_pairs_N.pdf
manuscript/figures/prospect_pairs_Cab.png: manuscript/figures/prospect_pairs_Cab.pdf
manuscript/figures/prospect_pairs_Car.png: manuscript/figures/prospect_pairs_Car.pdf
manuscript/figures/prospect_pairs_Cw.png: manuscript/figures/prospect_pairs_Cw.pdf
manuscript/figures/prospect_pairs_Cm.png: manuscript/figures/prospect_pairs_Cm.pdf

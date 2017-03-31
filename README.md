# tSDRG
Tree tensor network (TTN) strong disorder renormalisation group (SDRG) algorithm.

This is the code used in the paper "Self-assembling tensor networks and holography in disordered spin chains" by Andrew M. Goldsborough and Rudolf A. Römer, Phys. Rev. B 89, 214203. https://doi.org/10.1103/PhysRevB.89.214203 

The main scripts are:
- tSDRG_corr => just calculates correlation functions
- tSDRG_ee => calculates correlations and entanglement entropy by simple bipartition
- tSDRG_eeb => calculates correlations and entanglement entropy of a block in the centre of the chain

The supporting files in ./data_analysis perform averaging over multiple disorder realisations. These were used to create the graphs in the paper above.

For more information please refer to the paper above and the PhD thesis of Andrew M. Goldsborough: http://wrap.warwick.ac.uk/id/eprint/70003

If this code is used in full or part to produce further scientific publications please cite our original paper.

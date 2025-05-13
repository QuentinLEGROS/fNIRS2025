# fNIRS2025
Reproducible results for the paper 'Enhancing Near-Infrared Spectroscopy Analysis Using Ordinal Pattern Methods'

The scripts 'ex_perm_entropy.m' and 'ex_PRSA.m' generate respectively the figs 2 and 3 in the paper.

To reproduce the results presented in the paper, first download the nireject-benchmark data (github links https://www.spiedigitallibrary.org/journals/neurophotonics/volume-11/issue-4/045008/NiReject--toward-automated-bad-channel-detection-in-functional-near/10.1117/1.NPh.11.4.045008.full).
Then, run 'SVM.m' to train the model, and 'Comp_methods.m' to evaluate the performance of each approach on the dataset.

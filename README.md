# SO-SCL

Subject to license: "**GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf**"

C implementation of Soft-Output Successive Cancellation List Decoding, **SO-SCL**, in mex format so it can be run from MATLAB.

To compile, in MATLAB: mex -O SOSCL_mex.c

Three MATLAB sample simulations are included:

1) sim_product decodes product codes with SO-SCL as the component decoder. Rows and columns are processed in MATLAB and their decoding is *not parallelized*. Output is recored in results.
2) sim_blkSO_acc empirically evaluates the calibration, in the Brier Score sense, of the blockwise SO.
3) sim_BLER_UER evaluates the SO control of undetected error rate.
Output is recored in results.

The following should be cited in association with results from this code.

1) P. Yuan, K. R. Duffy & M. Médard. "Near-optimal generalized decoding of Polar-like codes.", IEEE ISIT, 2024. 
2) P. Yuan, K. R. Duffy & M. Médard. "Soft-output successive cancellation list decoding", IEEE Transactions on Information Theory, 71 (2), 1007–1017, 2025.

For further details on GRAND, see: https://www.granddecoder.mit.edu/

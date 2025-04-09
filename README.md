# SO-SCL

Subject to license: "**GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf**"

C implementation of Soft-Output Successive Cancellation List Decoding, **SO-SCL**, in mex format so it can be run from MATLAB.

To compile, in MATLAB: mex -O SOSCL_mex.c; mex -O polarTrans.c

Three MATLAB sample simulations are included:

1) sim_product decodes product codes with SO-SCL as the component decoder. Rows and columns are processed in MATLAB and their decoding is *not parallelized*. Output is recored in results.
2) sim_blkSO_acc empirically evaluates the calibration, in the Brier Score sense, of the blockwise SO.
3) sim_BLER_UER evaluates the SO control of undetected error rate.
Output is recored in results.

The following should be cited in association with results from this code.

1) P. Yuan, M. Medard, K. Galligan, K. R. Duffy. "Soft-output (SO) GRAND and Iterative Decoding to Outperform LDPC Codes". IEEE Trans. Wireless Commun., 2025.
2) P. Yuan, K. R. Duffy & M. Médard. "Near-optimal generalized decoding of Polar-like codes.", IEEE ISIT, 2024. 
3) P. Yuan, K. R. Duffy & M. Médard. "Soft-output successive cancellation list decoding", IEEE Transactions on Information Theory, 71 (2), 1007–1017, 2025.

In this implementation, a polar-like code is characterized by four parameters: n, k, frz, and dCons.

1) n denotes the code length, which must be a power of two. If a different code length is required, techniques such as puncturing or shortening may be employed. 
2) k represents the message length, also referred to as the code dimension.
3) frz is a logical array of length n that specifies which bits in the u vector are frozen. A value of true (or 1) in frz indicates a frozen bit, while a value of false (or 0) denotes an information bit. Thus, the array frz contains exactly k entries with a value of false (0).
4) dCons is a matrix with two columns, specifying the constraints for dynamic frozen bits. For example, the matrix
[9, 8; 11, 8; 17, 14; 17, 16] defines three dynamic frozen bits: u_9 = u_8, u_11 = u_8, and u_17 = u_14 + u_16.

The function "polarTrans" performs the non-bit-reversal polar transform, mapping a length-n bit vector u to a length-n bit vector c via the transformation c = u*Fn, where Fn is the Kronecker power of order log2(n) of the matrix [1, 0; 1, 1].

The function preencode_dpolar maps a length-k message vector to a length-n bit vector u, according to the specifications given by frz and dCons. This process is also referred to as the pre-transform.

For further details on GRAND, see: https://www.granddecoder.mit.edu/

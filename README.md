# CSL 2023 list decoding
Low-rate tail-biting serial list decoder that uses a dual trellis

Looking at Jacob's polar vs TBCC paper we can analytically plot an accurate union bound of the  FER curve  for an infinite list size (it is FER because there are no erasures) from the determination of the entire distance spectrum.  We can also see from Jacob's paper that the total failure rate (FER + erasures) gets worse as list size gets shorter.  Real systems cannot use an infinite list size, and so they need to choose the "right" finite list size.  Our goal is to write a paper that directly addresses this design need.

As of today, there is no ANALYTICAL way to determine performance FER or erasure or total failure rate (TFR) as a function of SNR.  However, we have a conjecture (See Dec. 7 10:03 a.m.) about how to do this.  The specific task here is to take Jacob's code, simulate using SLVD at several specific maximum list sizes (probably powers of 2) and then develop analytical bounds that match the simulations.

As a closing argument in the paper, I think we can propose even with infinite complexity that there is a best choice of list size for performance that is smaller than infinite, which is the choice that maximizes capacity.

Specific tasks:
1. Implement serial list Viterbi decoder for Jacob's (32,512) CRC-TBCC.
transmitting the all-zeros codeword and with no noise determine the distance to each not-necessarily-tail-biting codeword as a function of list size.
2. Use these for upper and lower bounds on TFR as we discussed.
3. Compare with SLVD simulations to see if our bound is accurate.
4. Note the list size at which we find another CRC-TBCC codeword.

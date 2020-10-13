# Multiple Sequence Alignment operator

##### Description

`dist_alignment` computes pairwise distances between aligned biological sequences.

##### Usage

Input projection|.
---|---
`row`        |  factor, sequence names/IDs
`col`        |  numeric, position
`y-axis`        |  numeric, value corresponding to amino or nucleic acid
`color`        |  factor, optional, letter

Properties|.
---|---
`sequence_type` | whether it is `dna`, `rna`, or `protein` sequences
`matrix` | the matrix distance to be used, can be "similarity" or "identity"
`gap` | with nucleotides, if set to 1, gaps will be counted in the identity measure

Output relations|.
---|---
`dist_to`        | numeric, sequence name
`dist`        | numeric, distance

##### Details

This operator compute a matrix of pairwise distances from aligned sequences using similarity (Fitch matrix, for protein sequences only) or identity matrix (for protein and DNA sequences). The resulting matrix contains the squared root of the pairwise distances. For example, if identity between 2 sequences is 80 the squared root of (1.0 - 0.8) i.e. 0.4472136.

##### See Also

https://github.com/tercen/readfasta_operator

https://github.com/tercen/msa_operator

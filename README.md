# PMNS-and-redundancy
This repository contains code to generate PMNS, for efficient modular arithmetic (with randomisation if desired). It also allows to perform some analysis such as redundancy, representation distributions. The implementations are done using SageMath library http://www.sagemath.org/.
<br />
<br />
Here are the descriptions of the files in this repository:
 - pmns_generator.sage: contains functions to generate efficient PMNS, given a prime number.
 - EXAMPLE__pmns_gen.sage: shows an example of PMNS generation. 
 - pmns_arith_ops.sage: contains functions for arithmetic and conversion operations in the PMNS. It also contains all the tools to study the redundancy in the PMNS, unique representation computations and equality check.
 - EXAMPLE__redund.sage: presents examples of arithmetic and conversion operations, representation computations and equality check.
 - pmns_rand_ops.sage: contains functions for randomised arithmetic and conversion operations.
 - EXAMPLE__rando.sage: shows how to perform randomised operations.

Code title: "Riemannian preconditioning for tensor completion"

Implementation of the preconditioned nonlinear conjugate gradient algorithm for Tucker completion
(c) 2015-2016 Hiroyuki Kasai <kasai@is.uec.ac.jp> and  Bamdev Mishra <b.mishra@ulg.ac.be>

This package contains a MATLAB implementation of the algorithms presented in the report.

H. Kasai and B. Mishra,
"Riemannian preconditioning for tensor completion",
Technical report, arXiv:1506.02159, 2015.

This implementation is due to 
Hiroyuki Kasai <kasai@is.uec.ac.jp> and Bamdev Mishra <b.mishra@ulg.ac.be>, 2015.


The implementation is a research prototype still in development and is provided AS IS. 
No warranties or guarantees of any kind are given. Do not distribute this
code or use it other than for your own research without permission of the authors.

Feedback is greatly appreciated.


Installation:
-------------------------

- Set current directory as your current folder in Matlab or put it in your Matlab path.
- Run "Install_mex.m". You do not need to do this step for subsequent usage. Also, note that there are some precompiled mex files.
- Run "Run_me_first.m" to add folders to the working path. This needs to be done at the starting of each session.
- To check that everything works, run "Test.m" at Matlab command prompt
  (You should see some plots at the end.)


Files:
------
- Proposed_algo_files/fixedrankfactory_tucker_preconditioned.m: The factory description that contains the geometric characterization of the Tucker manifold.
- Proposed_algo_files/fixedrank_tensor_completion.m: The wrapper function with a list of options.
- Proposed_algo_files/tucker2multiarray.m: An additional implementation that converts a tensor in Tucker format, stored a structure, to full format which is now stored in a multiarray format. It is not really required for our purpose.


Disclaimer:
-----------

All contents are written by Hiroyuki Kasai and Bamdev Mishra except,

- "Auxiliary_files/*" are written by Michael Steinlechner, provided in the Matlab code package geomCG.
- "Manopt" is from the website http://manopt.org.


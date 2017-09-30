This code contains a variety of methods (see CompleteTensor.m) to complete the missing entries of a data tensor. Some of these methods assume a certain structure in the tensor, i.e. that the pattern of missing entries is 'column-structured', i.e. that the second dimension (in our case, genes) are either completely observed or completely missing.

This code also evaluates these methods using cross-validation (see TensorCV.m).

All method parameters used for the PSB paper are defined in GetArgs.m.

A good starting point are the scripts such as run_final_model and run_benchmarking.

To run unit tests, just type:
>> test_all


T = GetTensor('small_g978');

T2 = NormSigs(T);

assert(isequaln(T,T2));

T = GetTensor('small_g50');
T2 = NormSigs(T);

assert(~isequaln(T,T2));

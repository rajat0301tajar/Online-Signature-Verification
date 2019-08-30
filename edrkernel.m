function y=edrkernel(X,Y)
    e = 112.6499;
    edr_val = editdist(X,Y,e);
    alpha = 0.5;
    beta = 2;
    y = exp(-alpha*power(edr_val,beta));
end
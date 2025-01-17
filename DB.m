function [res] = DB(x)

[bp,bn] = bimu_bernoulli(x);

res=(bp./x).*(1-bn);

end


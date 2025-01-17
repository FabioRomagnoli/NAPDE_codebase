function Jp=Jacob_patt(l)

dd=ones(l,1);
d1=ones(l-1,1);
A=diag(dd)+diag(d1,1)+diag(d1,-1);
D=diag(dd);

Jp=[A,D,D;A,A,D;A,D,A];







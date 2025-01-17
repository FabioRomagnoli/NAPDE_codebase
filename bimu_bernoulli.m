
function [bp,bn] = bimu_bernoulli(x)
  % Restituisce bp=B(x) e bn=B(-x)

  % Check input
  if nargin ~= 1
    error("bimu_bernoulli: wrong number of input parameters.");
  end

  xlim= 1e-2;
  ax  = abs(x);
  bp  = zeros(size(x));
  bn  = bp;
  
  block1  = find(~ax);
  block21 = find((ax>80)&x>0);
  block22 = find((ax>80)&x<0);
  block3  = find((ax<=80)&(ax>xlim));
  block4  = find((ax<=xlim)&(ax~=0));
  
  %  X=0
  bp(block1)=1.;
  bn(block1)=1.;
  
  % ASYMPTOTICS
  bp(block21)=0.;
  bn(block21)=x(block21);
  bp(block22)=-x(block22);
  bn(block22)=0.;
  
  % INTERMEDIATE VALUES
  bp(block3)=x(block3)./(exp(x(block3))-1);
  bn(block3)=x(block3)+bp(block3); % è il termine con phi_i-phi_i+1, da molt all'u_i+1
  
  % SMALL VALUES
  if(any(block4))
    jj=1;
    fp=1.*ones(size(block4));
    fn=fp;
    df=fp;
    segno=1.;
    while (norm(df,inf) > eps)
      jj=jj+1;
      segno=-segno;
      df=df.*x(block4)/jj;
      fp=fp+df;
      fn=fn+segno*df;
    end
    bp(block4)=1./fp;
    bn(block4)=1./fn;
  end
  
    end

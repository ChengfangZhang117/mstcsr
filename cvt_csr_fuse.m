function [imgf,time] = cvt_csr_fuse(I1,I2,n,D,lambda, flag)
tic;
I1=double(I1);
I2=double(I2);

% n=3;
is_real = 1;

finest = 1 ;
nbscales=n;
C1 = fdct_wrapping(I1, is_real, finest, nbscales);
C2 = fdct_wrapping(I2, is_real, finest, nbscales);
[M,N]=size(I1);
 C=C1;
%  se=[4 4 4 4 4;4 16 16 16 4;4 16 64 16 4;4 16 16 16 4;4 4 4 4 4]/256;
for scalelevel=2:n
    for orientationlevel=1:length(C1{scalelevel})
        E1=C1{scalelevel}{orientationlevel};
        E2=C2{scalelevel}{orientationlevel};
%         E1=conv2(E1,se,'same');
%         E2=conv2(E2,se,'same');
%         W=abs(E1)>abs(E2);
    um=3;
    A1 = ordfilt2(abs(es2(E1,floor(um/2))), um*um, ones(um));
  	A2 = ordfilt2(abs(es2(E2,floor(um/2))), um*um, ones(um));
    % second step
  	W = (conv2(double(A1 > A2), ones(um), 'valid')) > floor(um*um/2);
        C{scalelevel}{orientationlevel}=C1{scalelevel}{orientationlevel}.*W+C2{scalelevel}{orientationlevel}.*~W;
    end
end

opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 500;
opt.rho = 100*lambda + 0.5;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
[X1, optinf1] = cbpdn(D, C1{1}{1}, lambda, opt);
[X2, optinf2] = cbpdn(D, C2{1}{1}, lambda, opt);

%activity level measure
A1=sum(abs(X1),3);
A2=sum(abs(X2),3);

if flag == 1  
    r=9;  
else
    r=3; 
end

ker=ones(2*r+1,2*r+1)/((2*r+1)*(2*r+1));
AA1=imfilter(A1,ker);
AA2=imfilter(A2,ker);
decisionMap=AA1>AA2;

%base layer fusion
if flag == 1  
    C{1}{1}= C1{1}{1}.*decisionMap +  C2{1}{1}.*(1-decisionMap);
else
   C{1}{1}=( C1{1}{1} +  C2{1}{1})/2;
end

imgf = ifdct_wrapping(C, is_real, M, N);
time=toc;
imgf=uint8(imgf);
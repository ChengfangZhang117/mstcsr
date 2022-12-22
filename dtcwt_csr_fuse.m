function [imgf,time] = dtcwt_csr_fuse(I1,I2,N,D,lambda, flag)
tic;
I1=double(I1);
I2=double(I2);
% N=2;
[Y1,h1] = dtwavexfm2(I1,N,'legall','qshift_06');
[Y2,h2] = dtwavexfm2(I2,N,'legall','qshift_06');

% npd = 16;
% fltlmbd = 5;
% [s1_l, s1_h] = lowpass(s1, fltlmbd, npd);
% [s2_l, s2_h] = lowpass(s2, fltlmbd, npd);

% Compute representation
%-------------------------------------------------------------------------%
%                               low-pass fusion
%-------------------------------------------------------------------------%
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 500;
opt.rho = 100*lambda + 0.5;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
[X1, optinf1] = cbpdn(D, Y1, lambda, opt);
[X2, optinf2] = cbpdn(D, Y2, lambda, opt);

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
    Y=Y1.*decisionMap + Y2.*(1-decisionMap);
else
    Y=(Y1 + Y2)/2;
end
%-------------------------------------------------------------------------%
%                               high-pass fusion
%-------------------------------------------------------------------------%
%detail layer fusion
h=h1;
for i=1:N
%     W{i}=E1{i}>E2{i};
%     W{i}=medfilt2(W{i});
    for j=1:6
    %   W=abs(h1{i,1}(:,:,j))>abs(h2{i,1}(:,:,j));
        um=3;
        A1 = ordfilt2(abs(es2(h1{i,1}(:,:,j),floor(um/2))), um*um, ones(um));
        A2 = ordfilt2(abs(es2(h2{i,1}(:,:,j),floor(um/2))), um*um, ones(um));
        %second step
        W = (conv2(double(A1 > A2), ones(um), 'valid')) > floor(um*um/2);
        h{i,1}(:,:,j)=h1{i,1}(:,:,j).*W+h2{i,1}(:,:,j).*~W;
    %   h{2,1}(:,:,j)=h1{2,1}(:,:,j).*W1+h2{2,1}(:,:,j).*~W1;
    end
end
F= dtwaveifm2(Y,h,'legall','qshift_06');
time=toc;
imgf=uint8(F);
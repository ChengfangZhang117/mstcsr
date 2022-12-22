function [imgf,time] = dwt_csr_fuse(img1, img2, zt,D, lambda, flag)

s1=double(img1);
s2=double(img2);
[hei,wid]=size(s1);

tic;
% Highpass filter test image

X=cell(zt,4);                                                                                       
Y=cell(zt,4);                                               
Z=cell(zt,4);  
for i=1:zt
    [X{i,1},X{i,2},X{i,3},X{i,4}]=dwt2(s1,'db1','mode','per'); 
    s1=X{i,1};
    [Y{i,1},Y{i,2},Y{i,3},Y{i,4}]=dwt2(s2,'db1','mode','per'); 
    s2=Y{i,1};
end

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
[X1, optinf1] = cbpdn(D, X{zt,1}, lambda, opt);
[X2, optinf2] = cbpdn(D, Y{zt,1}, lambda, opt);

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
    Z{zt,1}=X{zt,1}.*decisionMap + Y{zt,1}.*(1-decisionMap);
else
    Z{zt,1}=(X{zt,1} + Y{zt,1})/2;
end
%-------------------------------------------------------------------------%
%                               high-pass fusion
%-------------------------------------------------------------------------%
%detail layer fusion
for i=zt:-1:1                           
    for j=2:4
        Z{i,j}=selc(X{i,j},Y{i,j},3);    
    end
end


for i=zt:-1:1
    if i>1
        Z{i-1,1}=idwt2(Z{i,1},Z{i,2},Z{i,3},Z{i,4},'db1','mode','per');
    else
        imgf = idwt2(Z{i,1},Z{i,2},Z{i,3},Z{i,4},'db1','mode','per');
    end
end
time = toc;
imgf=uint8(imgf);
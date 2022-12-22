function [imgf,time] = nsct_csr_fuse(img1, img2, D, lambda, flag)
tic;
s1=double(img1);
s2=double(img2);
[hei,wid]=size(s1);

tic;
pfilter = 'pyrexc' ;              
dfilter = 'vk' ; 

nlevels = [2 3 3 4] ;        
% pfilter = 'maxflat' ;              
% dfilter = 'dmaxflat7' ; 

coeffs_1 = nsctdec( s1, nlevels, dfilter, pfilter );
coeffs_2 = nsctdec( s2, nlevels, dfilter, pfilter );

%Bandpass subbands 
[m,n]=size(s1);
coeffs=coeffs_2;
for i=2:numel(nlevels)+1

    if nlevels(i-1)==0
        E1=abs(coeffs_1{i});
        E2=abs(coeffs_2{i});
        %             map=E1>E2;
         um=3;
    A1 = ordfilt2(abs(es2(E1,floor(um/2))), um*um, ones(um));
  	A2 = ordfilt2(abs(es2(E2,floor(um/2))), um*um, ones(um));
    % second step
  	map= (conv2(double(A1 > A2), ones(um), 'valid')) > floor(um*um/2);
        coeffs{i}(map)=coeffs_1{i}(map);
    else
        for j=1:(2^nlevels(i-1))
            E1=abs(coeffs_1{i}{j});
            E2=abs(coeffs_2{i}{j});
%             map=E1>E2;
    um=3;
    A1 = ordfilt2(abs(es2(E1,floor(um/2))), um*um, ones(um));
  	A2 = ordfilt2(abs(es2(E2,floor(um/2))), um*um, ones(um));
    % second step
  	map= (conv2(double(A1 > A2), ones(um), 'valid')) > floor(um*um/2);
            coeffs{i}{j}(map)=coeffs_1{i}{j}(map);
        end
    end
end

% Lowpass subband

opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 500;
opt.rho = 100*lambda + 0.5;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
[X1, optinf1] = cbpdn(D, coeffs_1{1}, lambda, opt);
[X2, optinf2] = cbpdn(D, coeffs_2{1}, lambda, opt);

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
    coeffs{1}=coeffs_1{1}.*decisionMap + coeffs_2{1}.*(1-decisionMap);
else
    coeffs{1}=(coeffs_1{1} + coeffs_2{1})/2;
end

imgf= nsctrec( coeffs, dfilter, pfilter ) ;
time=toc;
imgf=uint8(imgf);
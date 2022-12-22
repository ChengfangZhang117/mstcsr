close all;
clear all;
clc;
path(path,'test_images/')
addpath(genpath('toolbox'));
% Load dictionary
D_csr=load('dict.mat');  
method = {'dwt-csr-4','dtcwt-csr-4','cvt-csr-4','nsct-csr-4'};
params.method = method;
%key parameters
lambda=0.01; 
overlap = 6;                    
epsilon=0.1;
level=4;
flag=1;
tic;
test_image=mygetdirfiles('test_images');
test_imagecell=load_images(test_image);
index_a = 1:2:size(test_imagecell);
index_b = 2:2:size(test_imagecell);
A=cell(size(index_a,2),1);
B=cell(size(index_b,2),1);
k = 1;%%显示方法
for i = 1:numel(test_imagecell)/2
    
    f = test_image{index_a(i)};
    [p, n, x] = fileparts(f);
    params.p = p;
    params.n = n;
    params.x1 = x;

     if size(test_imagecell{index_a(i)},3)>1
    A{i}=double(rgb2gray(test_imagecell{index_a(i)}));
    B{i}=double(rgb2gray(test_imagecell{index_b(i)}));
    else
    A{i}=double(test_imagecell{index_a(i)});
    B{i}=double(test_imagecell{index_b(i)});
     end
     if rem(size(A{i},1),8) == 0 && rem(size(A{i},2),8) ==0
        block_num = 8;
     else 
        [maxC,minC]=GCDLCM(size(A{i},1),size(A{i},2));
        block_num = maxC;
     end

    
 %%%dwt-csr
        if k <= size(method,2)
    fprintf('第%d组图像使用%s方法结果：\n',i, method{k});
        else
        k = 1;
        end
    [y_F_dwt_csr{i},time_dwt_csr]=dwt_csr_fuse(A{i},B{i},level,D_csr.D,lambda,flag);
     k = k+1;
    %%%dtcwt-csr
        if k <= size(method,2)
    fprintf('第%d组图像使用%s方法结果：\n',i, method{k});
        else
        k = 1;
        end
    [y_F_dtcwt_csr{i},time_dtcwt_csr]=dtcwt_csr_fuse(A{i},B{i},level,D_csr.D,lambda,flag);
       k = k+1;
    %%%cvt-csr
        if k <= size(method,2)
    fprintf('第%d组图像使用%s方法结果：\n',i, method{k});
        else
        k = 1;
        end
    [y_F_cvt_csr{i},time_cvt_csr]=cvt_csr_fuse(A{i},B{i},level+1,D_csr.D,lambda,flag);
    k = k+1;
    %%%nsct-csr
        if k <= size(method,2)
    fprintf('第%d组图像使用%s方法结果：\n',i, method{k});
        else
        k = 1;
        end
    [y_F_nsct_csr{i},time_nsct_csr]=nsct_csr_fuse(A{i},B{i},D_csr.D,lambda,flag);
    k = 1;
%%%将各种方法所得到的融合图像放到同一个大的矩阵中
   result = cat(3, y_F_dwt_csr{i}, y_F_dtcwt_csr{i}, y_F_cvt_csr{i}, y_F_nsct_csr{i});
   conf.fusion_image{i} = {};
%%%%将融合图像写入到指定的文件夹里
   for j = 1:numel(method)
        conf.fusion_image{i}{j} = fullfile(p, 'results', [n sprintf('[%d-%s]', j, method{j}) x]);
        imwrite(uint8(result(:, :, j)), conf.fusion_image{i}{j});%%%%将各种方法的结果放到result文件夹中
    end
end
time = toc;

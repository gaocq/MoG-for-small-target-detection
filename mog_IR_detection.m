function mog_IR_detection(readPath, savePath, temporal_step, patch, param, lr_prior, mog_prior)
% This code is for our paper: Gao, Chenqiang, et al. "Infrared small-dim target detection based on Markov random field guided noise modeling." Pattern Recognition 76 (2018): 463-475.
% If you use this code for your work, please cite this paper. 
% @article{Gao2018Infrared,
%    author = {Gao, Chenqiang and Wang, Lan and Xiao, Yongxing and Zhao, Qian and Meng, Deyu},
%    title = {Infrared small-dim target detection based on Markov random field guided noise modeling},
%    journal = {Pattern Recognition},
%    volume = {76},
%    number = {Supplement C},
%    pages = {463-475},
%    month = {2018/04/01/},
%    year = {2018}
% }

%% usage
% Some modules are written by C++, and you should recompile them if you run
% it with errors, Like as follows:
% >> mex update_Z.cpp
% >> mex reconstructImage.cpp


%% Create folders for saving results
if ~exist([savePath '/background'])  
    mkdir([savePath '/background']);   
end
for i=1:param.mog_k
    if ~exist([savePath sprintf('/%d', i)])  
        mkdir([savePath sprintf('/%d', i)]);  
    end
end

%% Get all image file names, please make sure that the image file order is correct by this reading way.
filesdir = dir([readPath '/*.jpg']);

if isempty( filesdir )
    filesdir = dir( [readPath '/*.bmp'] );
end
if isempty( filesdir )
    filesdir = dir([readPath '/*.png']);
end
if isempty( filesdir )
    fprintf('\n There is no any image in the folder of %s', readPath);
    return;
end
% get all image file names into a cell array;
files = { filesdir.name };

%% begin to process images using mog based detection method.
% get the size of images
[m n] = size( imread([readPath '/' files{1}]) );
for t=1:temporal_step: length(files)-patch.length+1 
    %% read images
    I_orig = zeros(m, n, patch.length, 'uint8');
    for tt = 1 : patch.length
        fprintf([readPath '/' files{tt+t-1} '\n']);
        I_orig(:, :, tt) = imread([readPath '/' files{tt+t-1}]);
    end
    
    % construct the patch image
    [X, locations] = constructPatchImage(I_orig, patch);
    % preprocess patch image
    X = double(X);   X = X/255;
    meanx = mean2(X); X = X-meanx;
    [lr_model, mog_model] = mog_rpca_markov(X, param, lr_prior, mog_prior, patch.size, patch.size);
    L_patch = lr_model.U * lr_model.V';
    label = mog_model.label;
    
    % recover the low-rank matrix
    L_patch1 = L_patch - min(L_patch(:));
    backgroundImage = reconstructPatchImage(L_patch1, patch, [m n], locations);
    
    % save into disk
    for tt = 1:patch.length
        I = mat2gray(backgroundImage(:,:,tt));
        imwrite(I, [savePath '/background/' files{tt+t-1}]);
    end
    
    % recover the noise components
    E_patch = X - L_patch;
    E_patch = E_patch - min(E_patch(:));
    for i= 1:param.mog_k
        EE_patch = zeros(size(E_patch));
        index = find(label == i);
        EE_patch(index) = E_patch( index );
        E = reconstructPatchImage(EE_patch, patch, [m n], locations);

        % save component images into disk
        for tt = 1:patch.length
            imwrite(mat2gray(E(:,:,tt)), [savePath sprintf('/%d/', i) files{tt+t-1}]);
        end
    end
end

function [patch_image locations] = constructPatchImage(I_orig, param)
sp = param.step;
sz = param.size;
[m n c]= size(I_orig);
if c ~= param.length
    error('The length of images dose not equel to the param.length');
end
patch_image = [];
locations = [];
for i=1: sp: m-sz
    for j=1:sp:n-sz
        temp = I_orig(i:i+sz-1, j:j+sz-1,:);
        patch_image = [patch_image, temp(:)];
        locations = [locations;i,j];
    end
end

function image = reconstructPatchImage(patch_image, param, imageSize, locations)
%% This function uses the C++ implementation to speed up.
FrameM = imageSize(1);
FrameN = imageSize(2);
sz = param.size;
recordTotalLen = ceil(param.size/param.step)^2;
locations = locations-1;
[m n] = size(patch_image);
image = reconstructImage(patch_image, locations, sz, m, n, FrameM, FrameN, recordTotalLen);



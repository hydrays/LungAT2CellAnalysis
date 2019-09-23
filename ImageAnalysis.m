A = imread('yourfilename.tif');
%% If multi-channel is used in imaging, set the channel corresponding to AT2 marker
%A = A(:,:,2);

C = imgaussfilt(A, 2);
C(C>75) = 75;
C = single(C)./max(single(C(:)));
B = imgaussfilt(C, 2);
L2 = imsegkmeans(B,2,'NormalizeInput',true);
bmask = boundarymask(L2);
%% Use this if the K-means labeling for the bright region (AT2 cells) is 1
%L2 = double(L2) - 1;
%% Use this if the K-means labeling for the bright region (AT2 cells) is 2
L2 = -double(L2) + 2;
CC = bwconncomp(L2);

stats = regionprops(CC, 'Area','Centroid'); 
idx = find([stats.Area] > 120); 
BW2 = ismember(labelmatrix(CC), idx); 
n = CC.NumObjects;
island_size_list = zeros(1, n);
for i = 1:n
   island_size_list(i) = length(CC.PixelIdxList{i});
end
[n1, x1] = hist(log(island_size_list));

centroids = cat(1, stats.Centroid);
at2center = centroids(idx,:);

%% Cut-off radius (in pixels)
r = 500;
Mdl = KDTreeSearcher(at2center);
[Idx, D] = rangesearch(Mdl,at2center,r);
N = length(Idx);
V = zeros(N, 1);
for i=1:N
    i
    subnet = at2center(Idx{i}, :);
    k = 2;
    innerMdl = KDTreeSearcher(subnet);
    [innerIdx, innerD] = knnsearch(innerMdl, subnet, 'K', k);
    if (length(innerD)>=2) 
        V(i) = mean(innerD(:,2), 'omitnan');
    else
        V(i) = r;
    end
end

%% Normalize and set the upper bound
vv = (V - min(V))/(max(V)-min(V));
vv(vv>0.1) = 0.1;
scatter(at2center(:,1), at2center(:,2), 1, vv, 'filled');
axis equal;
%% print('-dpng', '-r600', 'output.png')

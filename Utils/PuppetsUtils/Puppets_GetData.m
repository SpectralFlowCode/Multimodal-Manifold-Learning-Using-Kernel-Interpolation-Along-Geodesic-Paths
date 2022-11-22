%% Load images
datafiles_cam1=fullfile('Datasets','PuppetsDataset','s1_%d.jpg');% Camera 1 image files name format, where %d is the serial number.
datafiles_cam2=fullfile('Datasets','PuppetsDataset','s2_%d.jpg');% Camera 2 image files name format, where %d is the serial number.

n = DataParams.NumberOfDatapoints;  
ResizeFac=DataParams.ResizeFac;

% A sample image, to compute get the size of the vectors.
str0=fullfile('Datasets','PuppetsDataset',sprintf('s1_%d.jpg',100001));
Itmp = im2double(imread(str0));
Itmp = imresize(Itmp,ResizeFac);
m = size(Itmp, 1)*size(Itmp,2)*size(Itmp,3);
s1 = zeros(n, m);
s2 = zeros(n, m);

% A sample cropped image, to compute get the size of the vectors.
str0=fullfile('Datasets','PuppetsDataset',sprintf('s1_%d.jpg',100001));
Itmp = im2double(imread(str0));
Itmp=Itmp(:,1+(size(Itmp,2)/2):end,:);
Itmp = imresize(Itmp,ResizeFac);
m_chopped = size(Itmp, 1)*size(Itmp,2)*size(Itmp,3);
s_yoda = zeros(n, m_chopped);
s_bulldog = zeros(n, m_chopped);
s_bunny= zeros(n, m_chopped);


disp('Loading images for the puppets toy example');
reverseStr = '';
for i=1:n
    %i
    percentDone = 100 * i / n;
    msg = sprintf('   percentage done: %3.1f', percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    idx = 100010 + i;         % Index of sample to be loaded.
    str1=fullfile('Datasets','PuppetsDataset',sprintf('s1_%d.jpg',idx));
    str2=fullfile('Datasets','PuppetsDataset',sprintf('s2_%d.jpg',idx));
    I1 = im2double(imread(str1));           % Load snapshot from Camera 1.
    I2 = im2double(imread(str2));           % Load snapshot from Camera 2.
    
    rI1 = imresize(I1,ResizeFac);            % Downsample snapshot 1
    rI2 = imresize(I2,ResizeFac);            % Downsample snapshot 2
    s1(i,:) = reshape(rI1, [1, m]);         % Reshape snapshot into a vector, and store it.
    s2(i,:) = reshape(rI2, [1, m]);         % Reshape snapshot into a vector, and store it.
    
    IYoda=I1(:,1:160,:);
    IBulldog=I1(:,161:end,:);
    IBunny=I2(:,161:end,:);
    IYoda = imresize(IYoda,ResizeFac);
    IBulldog = imresize(IBulldog,ResizeFac);
    IBunny = imresize(IBunny,ResizeFac);
    s_yoda(i,:) = reshape(IYoda, [1, m_chopped]);         % Reshape snapshot into a vector, and store it.
    s_bulldog(i,:) = reshape(IBulldog, [1, m_chopped]);         % Reshape snapshot into a vector, and store it.
    s_bunny(i,:) = reshape(IBunny, [1, m_chopped]);         % Reshape snapshot into a vector, and store it.

end
fprintf([reverseStr]);

%% Extract angles
d_yoda = pdist2(s_yoda,s_yoda,'euclidean');
d_bulldog = pdist2(s_bulldog,s_bulldog,'euclidean');
d_bunny = pdist2(s_bunny,s_bunny,'euclidean');

tmpD=d_yoda;
tmpA=exp(-tmpD/median(tmpD(:))/NormFac);%figure; imagesc(db(tmpA));
[Kyoda,~] =SingleIterationNorm(tmpA,1);
tmpD=d_bulldog;
tmpA=exp(-tmpD/median(tmpD(:))/NormFac);%figure; imagesc(db(tmpA));
[Kbulldog,~] =SingleIterationNorm(tmpA,1);
tmpD=d_bunny;
tmpA=exp(-tmpD/median(tmpD(:))/NormFac);%figure; imagesc(db(tmpA));
[Kbunny,~] =SingleIterationNorm(tmpA,1);

numNeighbors_refine=20;
d_refined = pdist2(Kbunny',Kbunny','euclidean');
[MapEmbdbunny, ~ , ~ , ~]=DiffusionMapsFromDistance(d_refined,1,numNeighbors_refine);
d_refined = pdist2(Kbulldog',Kbulldog','euclidean');
[MapEmbdbulldog, ~ , ~ , ~]=DiffusionMapsFromDistance(d_refined,1,numNeighbors_refine);
d_refined = pdist2(Kyoda',Kyoda','euclidean');
[MapEmbdyoda, ~ , ~ , ~]=DiffusionMapsFromDistance(d_refined,1,numNeighbors_refine);

% AngleYoda=atan(MapEmbdyoda(:,2)./(sqrt(MapEmbdyoda(:,3).^2+MapEmbdyoda(:,2).^2)));
% AngleBulldog=atan(MapEmbdbulldog(:,2)./(sqrt(MapEmbdbulldog(:,3).^2+MapEmbdbulldog(:,2).^2)));
% AngleBunny=atan(MapEmbdbunny(:,2)./(sqrt(MapEmbdbunny(:,3).^2+MapEmbdbunny(:,2).^2)));

AngleYoda=atan(MapEmbdyoda(:,2)./(sqrt(MapEmbdyoda(:,3).^2+MapEmbdyoda(:,2).^2)+MapEmbdyoda(:,3)));
AngleBulldog=atan(MapEmbdbulldog(:,2)./(sqrt(MapEmbdbulldog(:,3).^2+MapEmbdbulldog(:,2).^2)+MapEmbdbulldog(:,3)));
AngleBunny =atan(MapEmbdbunny(:,2)./(sqrt(MapEmbdbunny(:,3).^2+MapEmbdbunny(:,2).^2)+MapEmbdbunny(:,3)));

% figure(); scatter(MapEmbdyoda(:,2),MapEmbdyoda(:,3),10,AngleYoda)
% figure(); scatter(MapEmbdbunny(:,2),MapEmbdbunny(:,3),10,AngleBunny)
% figure(); scatter(MapEmbdbulldog(:,2),MapEmbdbulldog(:,3),10,AngleBulldog)


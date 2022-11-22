N=DataParams.NumberOfDatapoints;

%% load dataset
if not(exist('test')==1) || not(exist('training')==1)
    load('Datasets\mnist.mat')
end
X=reshape(training.images(:,:,1:N),[],N);
labels=1+training.labels(1:N);

%% View1 - Rotations
X1=zeros(784,N);
n1=[];
for ind=1:N
    n1=[n1,2*DataParams.MaxRotation*rand()];
    tmp=imrotate(reshape(X(:,ind),28,28),n1(end),'crop');
    X1(:,ind)=tmp(:);
end

%% View2 - Occludeding rectangle 
locs=[];
for l_ind=1:numel(labels)
    locs=[locs,randsample(find(training.labels==(labels(l_ind)-1)),1)];
end
X2=reshape(training.images(:,:,locs),[],N);

n2=[];
for ind=1:N
    loc=1+mod(ind,28-DataParams.RectSize-1);
    n2=[n2,loc];
    tmp=reshape(X2(:,ind),28,28);
    tmp(loc:(loc+DataParams.RectSize),...
        loc:(loc+DataParams.RectSize))=1;
    X2(:,ind)=tmp(:);
end

%% Show dataset (sanity)
if 0
    figure();
    for ind =1:4
        tmp=randi(N);
        subplot(4,2,2*(ind-1)+1);        imagesc(reshape(X1(:,tmp),28,28))
        subplot(4,2,2*ind);              imagesc(reshape(X2(:,tmp),28,28))
    end
%     SaveFig(gcf,OutputFolder,FigPreamble+"_Illustration");
end

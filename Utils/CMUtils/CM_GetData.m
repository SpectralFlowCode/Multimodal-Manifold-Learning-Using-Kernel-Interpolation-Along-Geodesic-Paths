%% Download and unzip condition monitoring dataset
if ~exist(fullfile(pwd,'Datasets'),'dir')
    mkdir(fullfile(pwd,'Datasets'))
end
if ~exist(fullfile(pwd,'Datasets','CMDataset'),'dir')
    disp('Downloading condition monitoring dataset')
    if ~exist('CMDataset.zip','file')
        outfilename = websave(fullfile('Datasets','CMDataset.zip'),DatasetsParams.DownloadLink);
        disp(sprintf('Finished downloading condition monitorin dataset to: %s',fullfile(pwd,'GasDataset.zip')))
    end
    disp('Unzipping condition monitoring dataset''s folder')
    unzip(fullfile('Datasets','CMDataset.zip'),fullfile('Datasets','CMDataset'))
    disp(sprintf('Finished unzipping condition monitoring dataset to: %s',fullfile('Datasets','CMDataset')))
end


%% Parse dataset and organize data into cell arrays
SensorsToLoad={'PS1','PS2','PS3','PS4','PS5','PS6','TS1','TS2','TS3','TS4','EPS1','FS1','FS2','VS1','SE','CP','CE'};

if ~exist(fullfile(pwd,'Datasets','CMDataset','ParsedData.mat'))
    disp('Parsing data');
    Basefolder=fullfile(pwd,'Datasets','CMDataset');
    Samples=containers.Map;
    reverseStr = '';
    for ind=1:numel(SensorsToLoad)
        percentDone = 100 * ind / numel(SensorsToLoad);
        msg = sprintf('   percent done: %3.1f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        Samples(SensorsToLoad{ind}) = dlmread(fullfile(Basefolder,SensorsToLoad{ind}+".txt"))';
    end
    fprintf([reverseStr]);
    
    tmp=dlmread(fullfile(Basefolder,"profile"+".txt"))';
    RelevantFlags=tmp(end,:);
    
    ConditionNames={'Cooler','Valve','Pump','Accumulator'};
    Conditions=containers.Map;
    for ind=1:numel(ConditionNames)
        Conditions(ConditionNames{ind})=tmp(ind,:);
    end
    RelevantInds=find(RelevantFlags==0);    
    save(fullfile(pwd,'Datasets','CMDataset','ParsedData.mat'),'Samples','ConditionNames','Conditions','RelevantInds');
    disp(sprintf('Finished parsing, parsed data was saved to: %s',fullfile(pwd,'Datasets','CMDataset','ParsedData.mat')));
else
    disp(sprintf('Loading parsed data from %s',fullfile(pwd,'Datasets','CMDataset','ParsedData.mat')))
    load(fullfile(pwd,'Datasets','CMDataset','ParsedData.mat'));
end




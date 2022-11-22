%% Download and unzip dataset
if ~exist(fullfile(pwd,'Datasets'),'dir')
    mkdir(fullfile(pwd,'Datasets'))
end
if ~exist(fullfile(pwd,'Datasets','GasDataset'),'dir')
    disp('Downloading E-Nose dataset')
    if ~exist('GasDataset.zip','file')
        outfilename = websave(fullfile('Datasets','GasDataset.zip'),DatasetsParams.DownloadLink);
        disp(sprintf('Finished downloading E-Nose dataset to: %s',fullfile(pwd,'GasDataset.zip')))
    end
    disp('Unzipping E-Nose dataset''s folder')
    unzip(fullfile('Datasets','GasDataset.zip'),fullfile('Datasets','GasDataset'));
    disp(sprintf('Finished unzipping E-Nose dataset to: %s',fullfile(pwd,'Datasets','GasDataset')))
end

%% Parse dataset and save as mat file
if ~exist(fullfile(pwd,'Datasets','GasDataset','ParsedData.mat'))
    disp('Parsing data');
    ParsedData=[];
    reverseStr = '';
    for BatchInd=1:10
        percentDone = 100 * BatchInd / 10;
        msg = sprintf('   percent done: %3.1f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        BatchPath=fullfile(pwd,'Datasets','GasDataset',"batch"+num2str(BatchInd)+".dat");
        opts = detectImportOptions(BatchPath);
        opts.DataLines=[1 Inf];opts.Delimiter = {' ',':',';'};
        data = readtable(BatchPath,opts,'ReadVariableNames',false);
        data=str2double(table2array(data(:,1:end)));
        ParsedData=[ParsedData;data];
    end
    fprintf([reverseStr]);
    save(fullfile(pwd,"Datasets","GasDataset","ParsedData.mat"),'ParsedData');
    disp(sprintf('Finished parsing, parsed data was saved to: %s',fullfile(pwd,"Datasets","GasDataset","ParsedData.mat")))
else
    disp(sprintf('Loading parsed data from %s',fullfile(pwd,"Datasets","GasDataset","ParsedData.mat")))
    load(fullfile(pwd,"Datasets","GasDataset","ParsedData.mat"));
end

%% organize data into cell arrays
SensorsToLoad={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'};
SamplesCell=cell(1,6);
for SampleCellInd=1:6
    SamplesCell{SampleCellInd}=containers.Map;
    RelevantRows=find(ParsedData(:,1)==SampleCellInd);
    RelevantData=ParsedData(RelevantRows,:);
    for ind=1:numel(SensorsToLoad)
        Offset=4+(ind-1)*16;
        SensorsInd=Offset + [0:2:14];
        SamplesCell{SampleCellInd}(SensorsToLoad{ind}) = RelevantData(:,SensorsInd)';
    end
end

GasNames={'Ethanol', 'Ethylene', 'Ammonia', 'Acetaldehyde', 'Acetone', 'Toluene'};
GasGT=containers.Map;
for ind=1:numel(GasNames)
    RelevantRows=find(ParsedData(:,1)==ind);
    GasGT(GasNames{ind})=ParsedData(RelevantRows,2);
end
RelevantInds=(1:1000);




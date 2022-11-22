function [ColumnStochasticK,K,NoisyData] = GetKernelFromNoisySensor(Samples,SensorToUse,PreProcess,RelevantInds,NormFac,noise,NoiseParams,PlotSamples)
Data=Samples(SensorToUse);
tmpData=Data(:);
t=linspace(0,1,numel(tmpData))';
switch noise             
    case 'saw'
        Apd=NoiseParams.Apd;
        fd=NoiseParams.fd;
        e=(tmpData-median(tmpData));e=sort(e);
        medstd=e(0.99*numel(e))-e(0.01*numel(e));
        n=t*(Apd*medstd);
        n=mod(n,Apd*medstd/fd)*fd;
    case 'sine'
        As=NoiseParams.As;
        f=NoiseParams.f;
        e=(tmpData-median(tmpData));e=sort(e);
        medstd=e(0.99*numel(e))-e(0.01*numel(e));
        n=(As*medstd)*sin(2*pi*(1:numel(tmpData))/(size(Data,1)*f))';
    case 'drift'
        Ad=NoiseParams.Ad;
        e=(tmpData-median(tmpData));e=sort(e);
        medstd=e(0.99*numel(e))-e(0.01*numel(e));
        n=t*(Ad*medstd);
    case 'none'
        n=t*0;
    case 'spikes'
        Percentage=NoiseParams.Percentage;
        Aspikes=NoiseParams.Aspikes;
        n=zeros(size(t));
        i=randsample(numel(t),round(Percentage*size(Data,2)));
        Amps=(Aspikes(2)-Aspikes(1))*ones(numel(i),1)+Aspikes(1);
        e=(tmpData-median(tmpData));e=sort(e);
        medstd=e(0.9*numel(e))-e(0.1*numel(e));
        n(i)=Amps*medstd;   
    case 'truncs'
        Percentage=NoiseParams.Percentage;
        n=zeros(size(t));
        i=randsample(numel(t),round(Percentage*size(Data,2)));
        e=median(tmpData);
        e=max(tmpData);
end
if strcmp(noise,'truncs')
    tmpDatanoised=tmpData;
    tmpDatanoised(i)=e;
    if PlotSamples
        figure(); plot(tmpDatanoised-mean(tmpDatanoised));
        title(SensorToUse);
    end
else
    tmpDatanoised=tmpData+n;
    if PlotSamples
        if  strcmp(noise,'WGN')
            fsize=2*size(Data,1);
            figure(); plot(tmpData(1:fsize)); hold on;plot(tmpData(1:fsize)+n(1:fsize));
            title(SensorToUse);     
        else
            figure();
            plot(tmpData-mean(tmpData)+n);
            title(SensorToUse); axis tight;
            figure();
            plot(tmpData-mean(tmpData)); hold on; plot(n)
            title(SensorToUse); axis tight;
        end        
    end
end
NoisyData=buffer(tmpDatanoised,size(Data,1));
NoisyData=NoisyData(:,RelevantInds);
[ColumnStochasticK,K] = GetKernelFromData(NoisyData,PreProcess,NormFac);
end


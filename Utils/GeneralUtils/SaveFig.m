function [] = SaveFig(f,folder,filename,same_size)
if nargin==3
set(gcf, 'Position', get(0, 'Screensize'));
end
% set(gcf, 'Position',[140,-140,1920,1080]);
saveas(f,fullfile(folder,filename+".jpg"));

end


function [] = modePartMatrix(mode,partMatrix,stateNames)
    colBlue = [56, 97, 163]/255;
    f = figure;
    f.Position = [50 50 600 600]; % Max for IEEE paper (687 x 844 pixels)
    box on;
    b = bar(abs(partMatrix(:,mode)),'FaceColor',colBlue);
    xticks([1:length(partMatrix)]);
    xticklabels(stateNames);
    xtickangle(90);
    ylabel('Magnitude of part. factor');

end
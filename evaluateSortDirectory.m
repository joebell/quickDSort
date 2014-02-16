function evaluateSortDirectory()

    baseName = '~/Desktop/Data/140213/';
    wildcard = 'RL140213_004_*.mat';
    plotSize = [4,5];
    
    colorList = ['r','g','b','m','c','y','k','w'];
    plottingFigure = figure();
    
    fileList = dir([baseName,wildcard]);

    plotN = 1;
    page = 1;
    for fileN = 1:4 %length(fileList)
        
        load([baseName,fileList(fileN).name]);
        figure(plottingFigure);
        subplot(plotSize(1),plotSize(2),plotN);
        
        % Evaluate the sort
        clusterNumbers = unique(data.spikeClusters);
        nClusters = length(clusterNumbers);
        for spikeN = 1:length(data.spikeSamples)
            snippetBounds = round(data.spikeSamples(spikeN) + [-1 1].*data.spikeWidth*data.sampleRate);
            if ((snippetBounds(1) >= 1) && (snippetBounds(2) <= length(data.dVdT)))
                snippet = data.dVdT(snippetBounds(1):snippetBounds(2));
                for clustNn = 1:nClusters
                    clustN = clusterNumbers(clustNn);
                    spikeScore(spikeN,clustN) = sum(abs(snippet - data.spikeAvg{clustN}))/sum(abs(data.spikeAvg{clustN}));
                end
            else
                for clustNn = 1:nClusters
                    clustN = clusterNumbers(clustNn);
                    spikeScore(spikeN,clustN) = 0;
                end
            end
        end
        
        for clustNn = 1:nClusters
        	clustN = clusterNumbers(clustNn);
            ix = find(data.spikeClusters == clustN);
            scatter(spikeScore(ix,1),spikeScore(ix,2),4,'filled','Marker','+',...
                    'MarkerEdgeColor',colorList(clustN),'MarkerFaceColor',colorList(clustN));
            hold on;
            maxLim = max([xlim(),ylim()]);
            xlim([0 maxLim]); ylim([0 maxLim]);
            line(xlim(),ylim(),'Color','k','LineStyle','--');
        end
        
                
        title({fileList(fileN).name;['  M:',num2str(data.sortMetric)]},'Interpreter', 'none');
        plotN = plotN + 1;
        pause(.1);
        
        
        % When the page is full, save it.
        if ((plotN > plotSize(1)*plotSize(2)) || (fileN == length(fileList)))
            plotN = 1;
            savePDF(['SpikeSortingP',num2str(page),'.pdf']);
            close all;
            if (fileN < length(fileList))
                plottingFigure = figure();
                page = page + 1;
            end
        end
        
    end
    
    
    
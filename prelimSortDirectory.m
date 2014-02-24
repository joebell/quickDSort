function prelimSortDirectory()

    baseName = '~/Desktop/Data/140213/';
    wildcard = 'RL140213_004_*.mat';
    % clusterStructure = {[0,1,3,4],[1,1,2]};
    clusterStructure = {[0,1,2]};
    plotSize = [4,5];
    
    plottingFigure = figure();

    fileList = dir([baseName,wildcard]);

    plotN = 1;
    page = 1;
    for fileN = 1:length(fileList)
    
    	load([baseName,fileList(fileN).name]);
    	if (true) %data.stimulus.stimNumber == 5

%        	unix(['cp ',baseName,fileList(fileN).name,' ','~/Desktop/ForRasters/',fileList(fileN).name]);
        
		    disp(fileList(fileN).name);
		    data = quickSort(data, clusterStructure);
		    [data, metric] = corrSortFast(data);
		    data = makeSpikeAvg(data);
            %[data, metric] = corrSortFast(data);
            %data = makeSpikeAvg(data);
		    
		    save([baseName,fileList(fileN).name],'data');
		    
            figure(plottingFigure);
		    subplot(plotSize(1),plotSize(2),plotN);
		    plotWaveforms(data); set(gca,'XTick',[],'YTick',[]);
		    title({fileList(fileN).name;['  M:',num2str(metric)]},'Interpreter', 'none');
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
    end

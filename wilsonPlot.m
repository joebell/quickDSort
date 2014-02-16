function wilsonPlot(scores,clusters)

    colorList = ['r','g','b','m','c','y','k'];

    nDim = size(scores,2);
    clusterList = unique(clusters);
    
    for dim1 = 1:nDim
        for dim2 = (dim1+1):nDim
            subplot(nDim-1,nDim-1,dim1+(dim2 - 2)*(nDim-1));
            for clustN = clusterList
                ix = find(clusters == clustN);
                scatter(scores(ix,dim1),scores(ix,dim2),4,'filled','Marker','+',...
                    'MarkerEdgeColor',colorList(clustN),'MarkerFaceColor',colorList(clustN));
                hold on;
                xlabel(['PC ',num2str(dim1)]);
                ylabel(['PC ',num2str(dim2)]);
                set(gca,'XTick',[],'YTick',[]);
            end
        end
    end
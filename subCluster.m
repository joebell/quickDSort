function T = subCluster(desc,labels)

    Nclust = length(labels);
	T = kmeans(desc,Nclust);
    % Sort the spikes by mean 1st descriptor
    for n = 1:Nclust
        ix = find(T == n);
        meanHeights(n) = mean(desc(ix,1));
    end
    
    [B,IX] = sort(abs(meanHeights),'descend');
    for N=1:Nclust
        n = IX(N);  ix = find(T == n);  Tp(ix) = labels(N);
    end
    T = Tp;
classdef MSPCA < handle
    properties
        accuracy;
        minTor;
    end
    
    properties(Access = private)
        MSName;
        MSMat;
        rawMSData;
        peakLocations;
        score;
    end
    
    properties(Dependent)
        nMS;
        nPeaks;
        data;
        pks;
    end
    
    properties(Dependent,Access=private)
        nMSCapacity;
        nMSLocCapacity;
    end
    
    methods
        function obj = MSPCA(accuracy,minTor)
            obj.accuracy = accuracy;
            obj.minTor = minTor;
            obj.MSName = {};
            obj.rawMSData = {};
            obj.MSMat = zeros(15,100);
            obj.peakLocations = zeros(1,100);
            obj.score = [];
        end
        
        function nMSCap = get.nMSCapacity(obj)
            nMSCap = size(obj.MSMat,1);
        end
        
        function nMSLocCap = get.nMSLocCapacity(obj)
            nMSLocCap = size(obj.MSMat,2);
        end
        
        function n_MS = get.nMS(obj)
            n_MS = length(obj.MSName);
        end
        
        function n_Peaks = get.nPeaks(obj)
            n_Peaks = sum(obj.peakLocations > 0);
        end
        
        function x = get.data(obj)
            x = obj.MSMat(1:obj.nMS,1:obj.nPeaks);
        end
        
        function pks = get.pks(obj)
            pks = obj.peakLocations(1:obj.nPeaks);
        end
        
        function addMS(obj,MSLoc,MSInts,name)
            try
                obj.rawMSData{obj.nMS+1} = [MSLoc,MSInts];
                obj.MSName{obj.nMS + 1} = name;
                [validInts,validLoc] = obj.findValidPeaks(MSInts,MSLoc);
                validLoc = round(validLoc/obj.accuracy) * obj.accuracy;
                for m = 1:1:length(validInts)
                    pkLoc = validLoc(m);
                    index = find(abs(obj.peakLocations - pkLoc)<obj.accuracy,1);
                    if isempty(index)
                        obj.setNewMSInts(obj.nPeaks+1,pkLoc,validInts(m));
                    else 
                        obj.setNewMSInts(index,pkLoc,validInts(m));
                    end
                end
            catch
                disp(strcat('Unable to add MS',32,name));
                obj.rawMSData(end) = [];
                obj.MSName(end) = [];
            end
        end
        
        function addMSFile(obj)
            [msCell,nameCell] = MatFile2MSs();       
            L = length(msCell);
            fprintf(1,'Add %d MS\n',L);
            for m = 1:1:L
                obj.addMS(msCell{m}(:,1),msCell{m}(:,2),nameCell{m});
            end          
        end
        
        function [raw,pksLoc,pksInts] = getMSById(obj,index)
            if index > obj.nMS
                return;
            else
                raw = obj.rawMSData{index};
                pksLoc = obj.peakLocations(1:obj.nPeaks);
                pksInts = obj.MSMat(index,1:obj.nPeaks);
            end
        end
        
        function [raw,pksLoc,pksInts] = getMSByName(obj,name)
            index = find(strcmp(obj.MSName,name));
            if isempty(index)
                fprintf(1,'No MS found named %s',name);
            else
                [raw,pksLoc,pksInts] = obj.getMSById(index);
            end
        end
        
        function sortMS(obj)
            [tmp,I] = sort(obj.pks);
            obj.peakLocations(1:obj.nPeaks) = tmp;
            tmp = obj.data(:,I);
            obj.MSMat(1:obj.nMS,1:obj.nPeaks) = tmp;
        end
        
        function [coeff,score,latent] = plotPCA(obj,isNormal,minInfo)
            obj.sortMS();
            tag = nameList2tags(obj.MSName);
            x = obj.MSMat(1:obj.nMS,1:obj.nPeaks);
            %%
            %x = x - repmat(mean(x),obj.nMS,1)./repmat(std(x),obj.nMS,1);
            if isNormal
                x = x./repmat(max(x,[],2),1,obj.nPeaks);
            end
            %%
            [coeff,score,latent] = pca(x);
            textPos = score + repmat(max(abs(score))*0.02,obj.nMS,1);
            figure;
            scatter(score(:,1),score(:,2),10,'filled');
            title(strcat('Total info:',32,num2str(sum(latent(1:2))/sum(latent))));
            box on;
            for m = 1:1:obj.nMS
                text(textPos(m,1),textPos(m,2),obj.MSName{m});
            end
            xlabel('PC1');
            ylabel('PC2');
            figure;
            scatter3(score(:,1),score(:,2),score(:,3),10,'filled');
            box on;
            title(strcat('Total info:',32,num2str(sum(latent(1:3))/sum(latent))));
            xlabel('PC1');
            ylabel('PC2');
            zlabel('PC3');
            for m = 1:1:obj.nMS
                text(textPos(m,1),textPos(m,2),textPos(m,3),obj.MSName{m});
            end
            
            figure;
            for m = 1:1:max(tag)
                scatter(score(tag==m,1),score(tag==m,2),10,'filled');
                hold on;
            end
            box on;
            xlabel('PC1');
            ylabel('PC2');
            title(strcat('Total info:',32,num2str(sum(latent(1:2))/sum(latent))));
            figure;
            for m = 1:1:max(tag)
                scatter3(score(tag==m,1),score(tag==m,2),score(tag==m,3),10,'filled');
                hold on;
            end
            box on;
            xlabel('PC1');
            ylabel('PC2');
            zlabel('PC3');
            title(strcat('Total info:',32,num2str(sum(latent(1:3))/sum(latent))));
            
            for m = 1:1:obj.nPeaks
                if(sum(latent(1:m))/sum(latent) > minInfo)
                    fprintf(1,'%d dimension for %.3f info\n',m,sum(latent(1:m))/sum(latent));
                    obj.score = score;
                    return;
                end
            end
               
        end
        
        function selectDim(obj,maxDim,selectDim,tryTime,GroupTag)
            if isempty(obj.score)
                disp('empty score!');
                return;
            end
            groupIndex = 1;
            for m = 1:1:tryTime
                dim = randsample(1:1:maxDim,selectDim);
                points = obj.score(:,dim);
                [~,gd] = groupDegree(points,GroupTag);
                if groupIndex > gd
                    fprintf(1,'GroupIndex: %.3f\n',gd);
                    disp(dim);
                    groupIndex = gd;
                end
            end         
        end
        
        function coefBar(obj,coeff,markTor)
            cm = [1,0,0;0,0,1];
            figure;
            hBar = bar(coeff(:,1:2));
            hBar(1).FaceColor = cm(1,:);
            hBar(2).FaceColor = cm(2,:);
            hBar(1).DisplayName = 'PC1';
            hBar(2).DisplayName = 'PC2';
            ids = cell(2,1);
            if markTor < 1
                for m = 1:1:2
                    ids{m} = find(abs(coeff(:,m))>max(abs(coeff(:,m)))*markTor);
                end
            elseif mod(markTor,1)==0
                for m = 1:1:2
                    [~,I] = sort(abs(coeff(:,m)),'descend');
                    ids{m} = I(1:markTor);
                end
            else
                disp('Invalid markTor!');
            end
            hold on;
            for m = 1:2
                pIndex = ids{m};
                L = length(pIndex);
                for n = 1:1:L
                    text(pIndex(n),coeff(pIndex(n),m),num2str(obj.pks(pIndex(n))),...
                        'FontSize',8,'Color',cm(m,:));
                end
            end
        end
        
        function plotMSByName(obj,name,hAxes)
            [raw,~,~] = obj.getMSByName(name);
            if ~exist('hAxes','var')
                figure; hAxes = axes;
            end
            plot(hAxes,raw(:,1),raw(:,2));
            title(name);
        end
        
    end
    
    methods(Access = private)
        function [nC,C] = getDangerCluster(obj,tag)
            L = length(tag);
            tmp = 0;
            C = {};
            nC = 0;
            for m = 1:1:L
                switch(2*tmp - tag(m))
                    case -1
                        %01
                        nC = nC + 1;
                        C{nC} = [m];
                    case 1
                        %11
                        C{nC}(end+1) = m;
                end
                tmp = tag(m);
            end
            
            for m = 1:1:nC
                C{m}(end+1) = C{m}(end) + 1;
            end
        end
        function [pksInts,pksLoc] = findValidPeaks(obj,MSInts,MSLoc)
            if any((MSLoc(2:end) - MSLoc(1:(end-1)))<0)
                [MSLoc,I] = sort(MSLoc);
                MSInts = MSInts(I);
            end
            pos = find((MSLoc(2:end) - MSLoc(1:(end-1)))==0);
            if ~isempty(pos)
                for m = 1:1:length(pos)
                    d = MSLoc(pos(m)+2) - MSLoc(pos(m)+1);
                    MSLoc(pos(m)+1) = MSLoc(pos(m)+1) + min(0.0001,d/2);
                end
            end
            [pksInts,pksLoc] = findpeaks(MSInts,MSLoc);
            filter = pksInts >= (max(pksInts) * obj.minTor);
            pksInts = pksInts(filter);
            pksLoc = pksLoc(filter);
            peakSub = abs(pksLoc(1:(end-1)) - pksLoc(2:end));
            if any(peakSub < obj.accuracy)
                [nC,danger_clusters] = obj.getDangerCluster(peakSub < obj.accuracy);
                for m = 1:1:nC
                    [~,I] = max(pksInts(danger_clusters{m}));
                    danger_clusters{m}(I) = [];
                end
                toDelet = [];
                for m = 1:1:nC
                    toDelet = [toDelet,danger_clusters{m}];
                end
                pksLoc(toDelet) = [];
                pksInts(toDelet) = [];
            end
        end
        function setNewMSInts(obj,index,loc,ints)
            if obj.nMS > obj.nMSCapacity
                [r,c] = size(obj.MSMat);
                newMat = zeros(r+5,c);
                newMat(1:r,1:c) = obj.MSMat;
                obj.MSMat = newMat;
            end
            
            if index > obj.nMSLocCapacity
                [r,c] = size(obj.MSMat);
                newCap = round(c*1.5);
                newMat = zeros(r,newCap);
                newPeakLocations = zeros(1,newCap);
                newMat(1:r,1:c) = obj.MSMat;
                newPeakLocations(1,1:c) = obj.peakLocations;
                obj.MSMat = newMat;
                obj.peakLocations = newPeakLocations;
                obj.peakLocations(1,c+1) = loc;
                obj.onPeakAdding(c+1,loc);
                return;
            else if index > obj.nPeaks
                    c = obj.nPeaks;
                    obj.peakLocations(1,c+1) = loc;
                    obj.onPeakAdding(c+1,loc);
                    return;
                end
            end
            
            obj.MSMat(obj.nMS,index) = ints;
        end
        
        function onPeakAdding(obj,index,loc)
            for m = 1:1:obj.nMS
                raw = obj.rawMSData{m};
                filter = abs(raw(:,1) - loc) < obj.accuracy;
                if any(filter)
                    obj.MSMat(m,index) = max(raw(filter,2));
                else
                    obj.MSMat(m,index) = 0;
                end
            end
        end
    end
    
end


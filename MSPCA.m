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
        
        function subMSLoc(obj,MSLoc)
            MSLoc = round(MSLoc/obj.accuracy) * obj.accuracy;
            index = find(abs(obj.peakLocations - MSLoc)<obj.accuracy,1);
            if ~isempty(index)              
                fprintf(1,'Peak Loc: %.3f has been substrated\n',obj.peakLocations(index));
                obj.peakLocations(index)=[];
                obj.MSMat(:,index)=[];
            end
        end
        
        function addMSFile(obj,isView)
            if nargin == 1
                isView = 0;
            end
            [msCell,nameCell] = MatFile2MSs();       
            L = length(msCell);
            fprintf(1,'Add %d MS\n',L);
            for m = 1:1:L
                obj.addMS(msCell{m}(:,1),msCell{m}(:,2),nameCell{m});
                if isView
                    figure('Position',[0,0,800,200]);
                    plot(msCell{m}(:,1),msCell{m}(:,2));
                    title(nameCell{m});
                end
            end
            choice = questdlg('Do you want to substrate background?','MSPCA','Yes','No','No');
            if strcmp(choice,'Yes')
                [msCell,~] = MatFile2MSs();
                [~,locSet] = obj.findValidPeaks(msCell{1}(:,2),msCell{1}(:,1));
                if length(msCell) > 1
                    locSet = round(locSet/obj.accuracy)*obj.accuracy;
                    choice = questdlg('More than 1 columns were detected,which operation do you want?','MSPCA','intersect','union','intersect');
                    if strcmp(choice,'intersect')
                        for m = 2:1:length(msCell)
                            [~,tmp] = obj.findValidPeaks(msCell{m}(:,2),msCell{m}(:,1));
                            tmp = round(tmp/obj.accuracy)*obj.accuracy;
                            locSet = intersect(tmp,locSet);
                        end
                    elseif strcmp(choice,'union')
                        for m = 1:1:length(msCell)
                            [~,tmp] = obj.findValidPeaks(msCell{m}(:,2),msCell{m}(:,1));
                            tmp = round(tmp/obj.accuracy)*obj.accuracy;
                            locSet = union(tmp,locSet);
                        end
                    else
                        return;
                    end
                end
                for m = locSet'
                    obj.subMSLoc(m);
                end
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
                fprintf(1,'No MS found named %s\n',name);
                raw = [];
                pksLoc = [];
                pksInts = [];
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
        
        function [names,tag,dic] = getTag(obj)
            names = unique(obj.MSName);
            [tag,dic] = nameList2tags(obj.MSName);
        end
        
        function [res] = pPCA(obj,isNormal)
            res = struct(); res.method = 'PCA';
            obj.sortMS();
%             tag = nameList2tags(obj.MSName);
            x = obj.MSMat(1:obj.nMS,1:obj.nPeaks);
            %%
            if isNormal
                x = x./repmat(max(x,[],2),1,obj.nPeaks);
            end
            %%
            [coe,sco,lat] = pca(x);
            res.score = sco;
            res.coeff = coe;
            res.latent = lat;
            obj.score = sco;
%             textPos = score + repmat(max(abs(score))*0.02,obj.nMS,1);
%             figure;
%             scatter(score(:,1),score(:,2),15,'filled');
%             title(strcat('Total info:',32,num2str(sum(latent(1:2))/sum(latent))));
%             box on;
%             for m = 1:1:obj.nMS
%                 text(textPos(m,1),textPos(m,2),obj.MSName{m});
%             end
%             xlabel(sprintf('PC1 (%.2f %%)',100*latent(1)/sum(latent)));
%             ylabel(sprintf('PC2 (%.2f %%)',100*latent(2)/sum(latent)));
%             figure;
%             scatter3(score(:,1),score(:,2),score(:,3),10,'filled');
%             box on;
%             title(strcat('Total info:',32,num2str(sum(latent(1:3))/sum(latent))));
%             xlabel('PC1');
%             ylabel('PC2');
%             zlabel('PC3');
%             for m = 1:1:obj.nMS
%                 text(textPos(m,1),textPos(m,2),textPos(m,3),obj.MSName{m});
%             end
%             
%             figure; ha = gca;
%             for m = 1:1:max(tag)
%                 scatter(score(tag==m,1),score(tag==m,2),10,'filled');
%                 %drawConfiInter(ha,score(tag==m,1:2));
%                 hold on;
%             end
%             %drawConfiInter(ha,score(tag==m,1:2));
%             box on;
%             xlabel(sprintf('PC1 (%.2f %%)',100*latent(1)/sum(latent)));
%             ylabel(sprintf('PC2 (%.2f %%)',100*latent(2)/sum(latent)));
%             title(strcat('Total info:',32,num2str(sum(latent(1:2))/sum(latent))));
%             figure;
%             for m = 1:1:max(tag)
%                 scatter3(score(tag==m,1),score(tag==m,2),score(tag==m,3),10,'filled');
%                 hold on;
%             end
%             box on;
%             xlabel('PC1');
%             ylabel('PC2');
%             zlabel('PC3');
%             title(strcat('Total info:',32,num2str(sum(latent(1:3))/sum(latent))));
%             
%             for m = 1:1:obj.nPeaks
%                 if(sum(latent(1:m))/sum(latent) > minInfo)
%                     fprintf(1,'%d dimension for %.3f info\n',m,sum(latent(1:m))/sum(latent));
%                     obj.score = score;
%                     return;
%                 end
%             end
               
        end
        
        function [res] = pTSNE(obj,isNormal,perp,dist)
            if ~exist('perp','var')
                perp = 10;
            end
            if ~exist('dist','var')
                dist = 'euclidean';
            end
            obj.sortMS();
%             tag = nameList2tags(obj.MSName);
            x = obj.MSMat(1:obj.nMS,1:obj.nPeaks);
            %%
            %x = x - repmat(mean(x),obj.nMS,1)./repmat(std(x),obj.nMS,1);
            if isNormal
                x = x./repmat(max(x,[],2),1,obj.nPeaks);
            end
            %%            
            Y = tsne(x,'Distance',dist,'Perplexity',perp);      
            res = struct();
            res.method = 'tSNE';
            res.score = Y;
            obj.score = Y;
        end
        
        function [dim,gd] = selectDim(obj,res,Dim,nSelectDim,tryTime,GroupTag)
            if isempty(obj.score)
                disp('empty score!');
                return;
            end
            dim = 1:nSelectDim;
            [~,gd] = groupDegree(res.score(:,dim),GroupTag);
            groupIndex = gd;
            for m = 1:1:tryTime
                tdim = randsample(1:1:Dim,nSelectDim);
                points = res.score(:,tdim);
                [~,tgd] = groupDegree(points,GroupTag);
                if groupIndex > tgd
                    fprintf(1,'GroupIndex: %.3f  ',tgd);
                    disp(tdim);
                    groupIndex = tgd;
                    gd = tgd;
                    dim = tdim;
                end
            end         
        end
        
        function coefBar(obj,res,markTor)
            if strcmp(res.method,'PCA')
                cm = [1,0,0;0,0,1];
                figure;
                hBar = bar(res.coeff(:,1:2));
                hBar(1).FaceColor = cm(1,:);
                hBar(2).FaceColor = cm(2,:);
                hBar(1).DisplayName = 'PC1';
                hBar(2).DisplayName = 'PC2';
                ids = cell(2,1);
                if markTor < 1
                    for m = 1:1:2
                        ids{m} = find(abs(res.coeff(:,m))>max(abs(res.coeff(:,m)))*markTor);
                    end
                elseif mod(markTor,1)==0
                    for m = 1:1:2
                        [~,I] = sort(abs(res.coeff(:,m)),'descend');
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
        end
        
        function coefScatter(obj,coeff,markTor)
            cm = [1,0,0;0,0,1];
            figure;
            scatter(coeff(:,1),coeff(:,2),10,'filled');
            xlabel('PC1'); ylabel('PC2');
            if markTor < 1
                ratio = abs(coeff(:,1:2))./repmat(max(abs(coeff(:,1:2))),size(coeff,1),1);
                I = find(max(ratio,[],2)>markTor);
            elseif mod(markTor,1)==0
                maxValue = max(abs(coeff(:,1:2)),[],2);
                [~,I] = sort(maxValue,'descend');
                I = I(1:markTor);
            else
                disp('Invalid markTor!');
            end
            hold on;
      
            L = length(I);
            for m = 1:1:L
                text(coeff(I(m),1),coeff(I(m),2),num2str(obj.pks(I(m))),...
                    'FontSize',8);
            end
            line([min(coeff(:,1))-0.5*range(coeff(:,1)),...
                max(coeff(:,1))+0.5*range(coeff(:,1))],[0,0],'Color','r','LineStyle','--');
            line([0,0],[min(coeff(:,2))-0.5*range(coeff(:,2)),...
                max(coeff(:,2))+0.5*range(coeff(:,1))],'Color','r','LineStyle','--');
            box on;
        end
        
        function plotMSByName(obj,name,hAxes)
            [raw,~,~] = obj.getMSByName(name);
            if ~exist('hAxes','var')
                figure; hAxes = axes;
            end
            plot(hAxes,raw(:,1),raw(:,2));
            title(name);
        end
        
        function scatterRes(obj,res,dim,CRegOpt)
            if ~exist('dim','var')
                dim = [1,2];
            end
            [~,tag,dic] = obj.getTag(); key = dic.keys;
            figure; ha = gca; ha.NextPlot = 'add';
            index = unique(tag); L = length(index);
            c = lines(L);
            for m = 1:L
                scatter(res.score(tag==index(m),dim(1)),res.score(tag==index(m),dim(2)),...
                    20,c(m,:),'filled','DisplayName',key{m});
                if exist('CRegOpt','var')
                    % ref: http://bbs.pinggu.org/thread-3216247-1-1.html
                    obj.drawConfReg(ha,res.score(tag==index(m),dim),c(m,:),...
                        CRegOpt.alpha,CRegOpt.dist);
                end
            end
            if strcmp(res.method,'PCA')
                xlabel(sprintf('PC1(%.2f%%)',res.latent(1)*100/sum(res.latent)));
                ylabel(sprintf('PC2(%.2f%%)',res.latent(2)*100/sum(res.latent)));
            end
            box on;
            
            figure; ha = gca; ha.NextPlot = 'add';
            for m = 1:L
                I = tag==index(m);
                scatter(res.score(I,dim(1)),res.score(I,dim(2)),...
                    20,c(m,:),'filled','DisplayName',key{m});
                text(res.score(I,dim(1)),res.score(I,dim(2)),obj.MSName(I),'Color',c(m,:));
            end
            box on;
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
        function drawConfReg(obj,ha,data,c,alpha,dist)
            if ~exist('alpha','var')
                alpha = 0.05;
            end
            if ~exist('dist','var')
                dist = 'norm';
            end
            ha.NextPlot = 'add';
            meanValue = mean(data);
            diffMat = data - repmat(meanValue,size(data,1),1);
            s = inv(cov(data));
            rd = sum(diffMat*s.*diffMat,2);
            if strcmp(dist,'norm')
                r = chi2inv(1-alpha,2);
            elseif strcmp(dist,'exp')
                r = prctile(rd,100*(1-alpha));
            else
                error('invalid distribution %s',dist);
            end
            % draw
            [V,D] = eig(s);
            aa = sqrt(r/D(1)); bb = sqrt(r/D(4));
            t = linspace(0,2*pi,100);
            xy = V*[aa*cos(t);bb*sin(t)]; xy = xy + repmat(mean(data)',1,100);
            plot(ha,xy(1,:),xy(2,:),'Color',c,'LineWidth',1);
        end
    end
    
end


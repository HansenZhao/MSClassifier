function [ tags,taggerDic ] = nameList2tags( nameCell,fitNum )
    L = length(nameCell);
    if ~exist('fitNum','var')
        fitNum = 1;
    end
    tags = ones(L,1);
    tmpTag = 1;
    tmpStr = nameCell{1}(1:fitNum);
    taggerDic = containers.Map();
    for m = 1:L
        newStr = nameCell{m}(1:fitNum);
        if ~taggerDic.isKey(newStr)
            taggerDic(newStr) = taggerDic.Count + 1;
        end
        tags(m) = taggerDic(newStr);
    end
    
end


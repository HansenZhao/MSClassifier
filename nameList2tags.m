function [ tags ] = nameList2tags( nameCell,fitNum )
    L = length(nameCell);
    if ~exist('fitNum','var')
        fitNum = 1;
    end
    tags = ones(L,1);
    tmpTag = 1;
    tmpStr = nameCell{1}(1:fitNum);
    for m = 2:L
        newStr = nameCell{m}(1:fitNum);
        if ~strcmp(tmpStr,newStr)
            tmpTag = tmpTag + 1;
            tmpStr = newStr;
        end
        tags(m) = tmpTag;
    end
    
end


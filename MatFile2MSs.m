function [ MScell,nameCell ] = MatFile2MSs()
    % convert MS1 file with many samples
    [fn,fp,index] = uigetfile('*.csv','please select MS file...');
    if index
        x = importdata(strcat(fp,fn));
        Mat = x.data;
        msName = x.textdata;
        [~,c] = size(Mat);
        msCount = (c+1)/3;
        if mod(msCount,1) == 0
            MScell = cell(msCount,1);
            nameCell = cell(msCount,1);
            for m = 1:1:msCount
                start = 3*(m-1) + 1;
                length = sum(Mat(:,start)>0);
                MScell{m} = Mat(1:length,start:(start+1));
                nameCell{m} = x.textdata{start};
            end
        else
            error('File Format Error!');
        end
    else
        error('No file selected!');
    end
end


function fdata = diagonal_filter(data, threshold)

struct = fliplr(eye(3));
fdata = normalize_for_filter(data, threshold);

background = imerode(fdata > 0, struct);

maxdata = ordfilt2(fdata, 9, true(3), struct);

fdata = fdata .* background;
image(10*fdata)
end

function data = normalize_for_filter(data, threshold)

% Remove outliers
% image(data)
% pause
% nzdata = data(data ~= 0);
% nzdata = reshape(nzdata, 1, numel(nzdata));
% q1 = quantile(nzdata, 0.25);
% q3 = quantile(nzdata, 0.75);
% iqr = q3-q1;
% data(data > q3 + 0.5*iqr) = q3;
% 
% image(data)
% pause

% Get rid of diagonal and bottom elements

offset = 8;
for i=1:size(data, 1)
    for j=i-offset:i+offset
        if(j < size(data, 1) && j > 0)
            data(i,j) = 0;
        end
    end
end
data(:,size(data,1)-5:end) = 0;
data(size(data,1)-5:end,:) = 0;

% Zscores
data(data < 0) = 0;
data = zscore(data, 0, 1);

% Make symmetric 
for i=1:size(data, 1)
    for j=i:size(data, 2)
        data(i,j) = data(j,i);
    end
end

% Threshold Zscores
data(data < threshold) = 0;

end

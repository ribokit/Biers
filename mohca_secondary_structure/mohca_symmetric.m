function res = mohca_symmetric(data)
res = data;
for i=1:size(data, 1)
    for j=i:size(data, 2)
        res(j,i) = res(i,j);
    end
end

end
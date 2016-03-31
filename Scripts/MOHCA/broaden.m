function res = broaden(data, n)
res = data;
struct = zeros(3);
struct(:, 2) = 1;
struct(2, :) = 1;
for i=1:n
    res = imdilate(res, struct);
end
end
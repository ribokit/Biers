function d_out = matrix_multiply_diagonal(d_in, x)

d_out = zeros(length(d_in),length(d_in));
mtx_tmp = zeros(length(d_in)+x-1,length(d_in)+x-1);
mtx_tmp((1+(x-1)/2):(length(mtx_tmp)-(x-1)/2),(1+(x-1)/2):(length(mtx_tmp)-(x-1)/2)) = d_in;

for i = 1:x
    d_out = d_out + mtx_tmp(i:(i+length(d_out)-1),i:(i+length(d_out)-1));    
end;
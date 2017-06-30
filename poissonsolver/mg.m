dims = 7

%V = [1 2 1 2 3 2 1 2 1   2 3 2 3 4 3 2 3 2   1 2 1 2 3 2 1 2 1]'
V = dims*ones(dims^3, 1);
%{
if (mod(size(V)(1), 2) == 0)
  V = V(1:end-1)
endif
%}

A = zeros(floor((size(V)(1))/27),size(V)(1));
B = zeros(size(V)(1), floor((size(V)(1))/27));

count = 1;
%{
for i = 2:2:dims
  for j = 2:2:dims
    for k = 2:2:dims
%}
    %restriction
    
const = dims^2+dims+1;
for i = 2:2:dims
  for j = 2:2:dims
    for k = 2:2:dims
      ind = i*dims*dims + j * dims + k + 1 - const;
      A(count, ind-dims^2-dims-1) = 1/64;
      A(count, ind-dims^2-dims) = 1/32;
      A(count, ind-dims^2-dims+1) = 1/64;
      A(count, ind-dims^2-1) = 1/32;
      A(count, ind-dims^2) = 1/16;
      A(count, ind-dims^2+1) = 1/32;
      A(count, ind-dims^2+dims-1) = 1/64;
      A(count, ind-dims^2+dims) = 1/32;
      A(count, ind-dims^2+dims+1) = 1/64;
      A(count, ind-dims-1) = 1/32;
      A(count, ind-dims) = 1/16;
      A(count, ind-dims+1) = 1/32;
      A(count, ind-1) = 1/16;
      A(count, ind) = 1/8;
      A(count, ind+1) = 1/16;
      A(count, ind+dims-1) = 1/32;
      A(count, ind+dims) = 1/16;
      A(count, ind+dims+1) = 1/32;
      A(count, ind+dims^2-dims-1) = 1/64;
      A(count, ind+dims^2-dims) = 1/32;
      A(count, ind+dims^2-dims+1) = 1/64;
      A(count, ind+dims^2-1) = 1/32;
      A(count, ind+dims^2) = 1/16;
      A(count, ind+dims^2+1) = 1/32;
      A(count, ind+dims^2+dims-1) = 1/64;
      A(count, ind+dims^2+dims) = 1/32;
      A(count, ind+dims^2+dims+1) = 1/64;
      %prolongation
      B(i, count) = 1/8;
      B(i+1, count) = 1/4;
      B(i+2, count) = 1/8;
      B(i+3, count) = 1/4;
      B(i+4, count) = 1/2;
      B(i+5, count) = 1/4;
      B(i+6, count) = 1/8;
      B(i+7, count) = 1/4;
      B(i+8, count) = 1/8;
      B(i+9, count) = 1/4;
      B(i+10, count) = 1/2;
      B(i+11, count) = 1/4;
      B(i+12, count) = 1/2;
      B(i+13, count) = 1;
      B(i+14, count) = 1/2;
      B(i+15, count) = 1/4;
      B(i+16, count) = 1/2;
      B(i+17, count) = 1/4;
      B(i+18, count) = 1/8;
      B(i+19, count) = 1/4;
      B(i+20, count) = 1/8;
      B(i+21, count) = 1/4;
      B(i+22, count) = 1/2;
      B(i+23, count) = 1/4;
      B(i+24, count) = 1/8;
      B(i+25, count) = 1/4;
      B(i+26, count) = 1/8;
      count = count + 1;
    endfor
  endfor
endfor

Vold = V;
V = A*V
V = B*V
[Vold V]
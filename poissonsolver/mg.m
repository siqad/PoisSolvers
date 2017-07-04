dims = 5

%V = [1 2 1 2 3 2 1 2 1   2 3 2 3 4 3 2 3 2   1 2 1 2 3 2 1 2 1]'
V = ones(dims^3, 1);
for i = 1:size(V)(1)
  if i > size(V)(1)/2
    V(i) = 2*V(i);
  endif
endfor
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
      count = count + 1;
    endfor
  endfor
endfor
B = (2^3)*A';

%Normalization
for i = 1:size(B)(1)
  B(i,:) = B(i,:) / sum(B(i,:));
endfor

Vold = V;
V = A*V
V = B*V
[Vold V]
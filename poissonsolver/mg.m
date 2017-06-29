%V = [1 2 1 2 3 2 1 2 1   2 3 2 3 4 3 2 3 2   1 2 1 2 3 2 1 2 1]'
V = [4 4 4 4 4 4 4 4 4   4 4 4 4 4 4 4 4 4   4 4 4 4 4 4 4 4 4   4 4 4 4 4 4 4 4 4   4 4 4 4 4 4 4 4 4]'

%{
if (mod(size(V)(1), 2) == 0)
  V = V(1:end-1)
endif
%}

A = zeros(floor((size(V)(1))/27),size(V)(1));
B = zeros(size(V)(1), floor((size(V)(1))/27));
i = 1;
count = 1;
while i < size(V)(1)
  %restriction
  A(count, i) = 1/64;
  A(count, i+1) = 1/32;
  A(count, i+2) = 1/64;
  A(count, i+3) = 1/32;
  A(count, i+4) = 1/16;
  A(count, i+5) = 1/32;
  A(count, i+6) = 1/64;
  A(count, i+7) = 1/32;
  A(count, i+8) = 1/64;
  A(count, i+9) = 1/32;
  A(count, i+10) = 1/16;
  A(count, i+11) = 1/32;
  A(count, i+12) = 1/16;
  A(count, i+13) = 1/8;
  A(count, i+14) = 1/16;
  A(count, i+15) = 1/32;
  A(count, i+16) = 1/16;
  A(count, i+17) = 1/32;
  A(count, i+18) = 1/64;
  A(count, i+19) = 1/32;
  A(count, i+20) = 1/64;
  A(count, i+21) = 1/32;
  A(count, i+22) = 1/16;
  A(count, i+23) = 1/32;
  A(count, i+24) = 1/64;
  A(count, i+25) = 1/32;
  A(count, i+26) = 1/64;
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
  i = i + 27;
  count = count + 1;
endwhile

Vold = V;
V = A*V
V = B*V
[Vold V]
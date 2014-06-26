#! /apps/OCTAVE/3.6.4/bin/octave -qf
#

arg_list = argv();
L = spconvert(load(arg_list{1}));
#L0 = spconvert(load(arg_list{2}));
load(arg_list{2});
[m,n] = size(L0);
L = L(1:m, 1:n);
R = norm(L - L0, inf)/ norm(L0, inf);
disp(R);

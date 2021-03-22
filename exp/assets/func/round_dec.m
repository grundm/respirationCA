function Y = round_dec(X,N)
% Y = round_dec(X,N) rounds X to nearest decimal digit position N
%
% Author:           Martin Grund
% Last update:      November 19, 2015

%%
Y = round(X*10^N)/10^N;
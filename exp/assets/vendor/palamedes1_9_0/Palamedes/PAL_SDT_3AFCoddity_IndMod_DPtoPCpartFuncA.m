%
% PAL_SDT_3AFCoddity_IndMod_DPtoPCpartFuncA
%   
%   Internal function
%
% Introduced: Palamedes version 1.6.0 (FK)

function Y = PAL_SDT_3AFCoddity_IndMod_DPtoPCpartFuncA(X,dP)

N = PAL_pdfNormal(X,0,1);
P = PAL_ZtoP(X+dP);
Y = N.*P.^2;
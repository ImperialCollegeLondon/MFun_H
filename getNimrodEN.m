function [E,N,EE,NN] = getNimrodEN()
%
% get EN location for NIMROD Raw Radar
% 5min
% 1725*2175

RD_info = 'K:\UK_Radar_Matlab\2013_12_29.mat';

load(RD_info,'DAT');

aaa = DAT.d201312290000;

dE = double(aaa.rl_gen_hd(6));
dN = double(aaa.rl_gen_hd(4));
E = double(aaa.rl_gen_hd(5));
N = double(aaa.rl_gen_hd(3));

[EE,NN] = meshgrid(E:dE:E+dE*(1725-1), N:-dN:N-dN*(2175-1));

end
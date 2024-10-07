function [H0, W, Grid_Sign] =  initialize_wave( param )
% function [H0, W, Grid_Sign] =  initialize_wave( param )
%
% This function return the wave height coefficients H0 and W for the
% parameters given in input. These coefficients are constants for a given
% set of input parameters.
% Third output parameter is optional (easy to recalculate anyway)

rng(param.rng);  %// setting seed for random numbers

gridSize = param.meshsize * [1 1] ;

meshLim = pi * param.meshsize / param.patchsize ;
N = linspace(-meshLim , meshLim , param.meshsize ) ;
M = linspace(-meshLim , meshLim , param.meshsize ) ;
[Kx,Ky] = meshgrid(N,M) ;

K = sqrt(Kx.^2 + Ky.^2);    %// ||K||
W = sqrt(K .* param.g);     %// deep water frequencies (empirical parameter)

[windx , windy] = pol2cart( deg2rad(param.winddir) , 1) ;

P = phillips(Kx, Ky, [windx , windy], param.windSpeed, param.A, param.g) ;
H0 = 1/sqrt(2) .* (randn(gridSize) + 1i .* randn(gridSize)) .* sqrt(P); % height field at time t = 0

if nargout == 3
    Grid_Sign = signGrid( param.meshsize ) ;
end
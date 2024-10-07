function Z = calc_wave( H0,W,time,Grid_Sign )
% Z = calc_wave( H0,W,time,Grid_Sign )
%
% This function calculate the wave height based on the wave coefficients H0
% and W, for a given "time". Default time=0 if not supplied.
% Fourth output parameter is optional (easy to recalculate anyway)

    % recalculate the grid sign if not supplied in input
    if nargin < 4
        Grid_Sign = signGrid( param.meshsize ) ;
    end
    % Assign time=0 if not specified in input
    if nargin < 3 ; time = 0 ; end

    wt = exp(1i .* W .* time ) ;
    Ht = H0 .* wt + conj(rot90(H0,2)) .* conj(wt) ;  
    Z = real( ifft2(Ht) .* Grid_Sign ) ;
end
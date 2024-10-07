function P = phillips(Kx, Ky, windDir, windSpeed, A, g)
%// The function now accept scalar, vector or full 2D grid matrix as input
    K_sq = Kx.^2 + Ky.^2;
    L = windSpeed.^2 ./ g;
    k_norm = sqrt(K_sq) ;
    WK = Kx./k_norm * windDir(1) + Ky./k_norm * windDir(2);
    P = A ./ K_sq.^2 .* exp(-1.0 ./ (K_sq * L^2)) .* WK.^2 ;
    P( K_sq==0 | WK<0 ) = 0 ;
end
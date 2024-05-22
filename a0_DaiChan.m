function [hn_LT] = a0_DaiChan(omega_L, omega_H, N)
    n = -(N-1)/2:(N-1)/2;
    hn_LT = zeros(1,length(n));
    for i = 1:length(n)
         if n(i) == 0
            hn_LT(i) = 1 - (omega_H - omega_L)/pi;
         else
            hn_LT(i) = -(sin(omega_H*n(i)) - sin(omega_L*n(i)))/(pi*n(i));
         end
    end
end
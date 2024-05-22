function [hn_LT] = a0_ThapQua(omega_c, N)
    n = -(N-1)/2:(N-1)/2;
    hn_LT = zeros(1,length(n));
    for i = 1:length(n)
         if n(i) == 0
            hn_LT(i) = omega_c/pi;
         else
            hn_LT(i) = sin(omega_c*n(i))/(pi*n(i));
         end
    end
end
function [hn_LT] = a0_DaDai(omega1_L, omega1_H, omega2_L, omega2_H, N)
    n = -(N-1)/2:(N-1)/2;
    hn_LT = zeros(1,length(n));
    for i = 1:length(n)
         if n(i) == 0
            hn_LT(i) = (omega1_H - omega1_L + omega2_H - omega2_L)/pi;
         else
            hn_LT(i) = (sin(omega1_H*n(i)) - sin(omega1_L*n(i)) + ...
                        sin(omega2_H*n(i)) - sin(omega2_L*n(i)))/(pi*n(i));
         end
    end
end
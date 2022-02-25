% generate the lowpass filter used in the DRNL implementation
% Morten Løve Jepsen, 2005

function [b,a]=lowpassDRNL(f0,fs);

    O=pi*f0/fs;
    
    C=1/(1+2^0.5*cot(O)+(cot(O))^2);
    D=2*C*(1-(cot(O))^2);
    E=C*(1-2^0.5*cot(O)+(cot(O))^2);

    b=[C, 2*C, C];
    a=[1, D, E];
    
 





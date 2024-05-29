function alpha_f=wateratten(freq)
% alpha_f: attenuation in water column [pos~db/m/kHz, neg~dB/wavelength]
% frequency in Hz
nfreq=length(freq); 
alpha_f=zeros(1,nfreq);
for ifreq=1:nfreq
    ff= freq(ifreq)/1000;  %ff -> frequency in kHz
    A=3.3e-03;
    B= 0.11*ff^2;
    C=1.0 + ff^2 ;
    D=44*ff^2 ;
    E=4100 + ff^2 ;
    F=(3.0e-04)*ff^2 ;
    %     alphaprime = A + B/C + D/E + F ;
    %     alphaprime=alphaprime/1000 ;
    %     alpha=alphaprime/8.686 ;
    %     alpha_f=8.686 * alpha/ff ;
    %     svp_in.walphs = alpha_f;
    alpha=(A + B/C + D/E + F )/8686 ;
    alpha_f(ifreq)=8.686 * alpha/ff ; %water attenuation in db/m/kHz
end

    

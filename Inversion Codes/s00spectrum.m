function [p,S]= s00spectrum(f,S_0)



    
  
    
S=S_0-10*log10(f.^3.594)+10*log10((1+(f/340).^2).^0.917); %Spectral Level in dB 



  % ind1995T2005=find( (f < 307) & ( f>=300) ); 
   
   
   %   S(ind1995T2005) = (66.6)*ones(size(ind1995T2005));

%p=10.^(S/20); % GOOD MOVE: removing uniform, random phase
randomPhz=rand(size(f))*2*pi-pi; %uniformly random phase from -pi to pi
p=10.^(S/20); % GOOD MOVE: removing uniform, random phase
 %p=10.^(S/20).*exp(1i*randomPhz); % GOOD MOVE: removing uniform, random phase


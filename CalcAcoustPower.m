function [Ea] = CalcAcoustPower(type,datavec, timevec, Specvar)
% calculate the acoustic power for a long range
% type  1: Method for local data with Waveform (filtered)
%       2: Method for remote data with single value (dB) with attenuation 
%            or '2b' single value without attenuation factor (Geometric spreading)
%               '2c' vector with attenuation 
%               '2d' vector without attenuation factor (geometric
%               spreading) 
% datavec = data type detecds on type
% timevec = vector of time stamps for the data 
% Specvar = structure witth the specific settings for the call 
%           ALL:     Specvar.Fac = factor 1,2,4
%           Type 1:  Specvar.R: distance from sensor to source in m, 
%                    Specvar.rho_air: (1 kg/m^3)air density,
%                    Specvar.cel: speed of sound in air m/sec, 
%                    Specvar.tau: time window (s),
%                    Specvar.fe: samplerate (1/sec)
%           Type 2: all of settings for type 1 and add 
%                    Specvar.Att for attenuation in dB

% Written by: Anna Perttu 
% Updated 07/03/2020 

isGlobal = islaSetGlobal;
switch type
    case 1
        Omega = Specvar.Fac*pi*Specvar.R^2; % geometric attenuation factor
        IntFac = 1/(Specvar.rho_air*Specvar.cel*Specvar.tau);
        
        % loop over tau
        t = timevec; % time vector
        ti= t(1);        %% initial and final times
        tf= max(t );
        
        x = datavec; % data vector (already detrended and filtered) 
        
        % 		everything should be integers from ti  to tf
        begwind=ti:Specvar.tau:tf;	%	Beginning of each window, in "step" increments.
             
        for j = 1:length(begwind)
            start = min(find(t  >= begwind(j))) ;  %% find starting index for window j
            endtime = max(find(t  <= begwind(j) + Specvar.tau)) ;  %% ends at max value before next window
            xint = x(start:endtime);  % select window of data from start to endtime
            data2 = xint.^2;
            Ea(j) = Omega*IntFac*(1/Specvar.fe).*trapz(data2); % acoustic power in Watts
        end
        
    case 2 % single value
       
        Omega = Specvar.Fac*pi*1000^2; % 1 km from vent! 
        IntFac = 1/(Specvar.rho_air*Specvar.cel);
        DelP = isGlobal.Aref*10^((datavec+Specvar.Att)/20);
        Ea = Omega*IntFac*DelP.^2;
    
    case '2b' % without attenuation factor 
        Omega = Specvar.Fac*pi*Specvar.R^2; % geometric attenuation factor
        IntFac = 1/(Specvar.rho_air*Specvar.cel);
        DelP = isGlobal.Aref*10^((datavec)/20);
        Ea = Omega*IntFac*DelP.^2;
    
    case '2c' % vector 
        
        Omega = Specvar.Fac*pi*1000^2; % 1 km from vent! 
        IntFac = 1/(Specvar.rho_air*Specvar.cel);
        
        t = timevec; % time vector
        ti= t(1);        %% initial and final times
        tf= max(t );     % 		everything should be integers from ti  to tf
        begwind=ti:Specvar.tau:tf;	%	Beginning of each window, in "step" increments.
        %
        for j = 1:length(begwind)
            start = min(find(t  >= begwind(j))) ;  %% find starting index for window j
            endtime = max(find(t  <= begwind(j) + Specvar.tau)) ;  %% ends at max value before next window
            xint = mean(datavec(start:endtime)); % mean in the time window changed 190218
            DelP = isGlobal.Aref*10.^((xint+Specvar.Att)./20);
            Ea(j,:) = Omega*IntFac*DelP.^2;
        end
        Ea(:,2) = (begwind);
    
    case '2d' % vector w/out calculated attenuation 
        
        Omega = Specvar.Fac*pi*Specvar.R^2; % 1 km from vent! 
        IntFac = 1/(Specvar.rho_air*Specvar.cel);
        
        t = timevec; % time vector
        ti= t(1);        %% initial and final times
        tf= max(t );     % 		everything should be integers from ti  to tf
        begwind=ti:Specvar.tau:tf;	%	Beginning of each window, in "step" increments.
        %
        for j = 1:length(begwind)
            start = min(find(t  >= begwind(j))) ;  %% find starting index for window j
            endtime = max(find(t  <= begwind(j) + Specvar.tau)) ;  %% ends at max value before next window
            xint = datavec(start:endtime);
            DelP = isGlobal.Aref*10.^((mean(xint))./20);
            Ea(j) = Omega*IntFac*DelP.^2;
        end    
        
end


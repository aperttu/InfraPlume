function [V_gas] = calcVelRadiationPat(Ea, Specvar);
% caclculate the vent gas exit velocity based on the acoustic power at the vent
% (in Watts) 
%
% Ea = acoustic power in Watts 
% Specvar = structure witth the specific settings for the call 
%                    Specvar.rho_air: (1 kg/m^3)air density,
%                    Specvar.cel: speed of sound in air m/sec, 
%                    Specvar.R_cond: conduit radius in m

% written by Anna Perttu
% updated 07/03/2020

for type=[1,2]
    if type == 1
        n = 4; m = 1; K = 1 ;	f = 4; % K=1(calcul Sylvie)	or K = 1/4 ; 4 is a question 8/21	%	Monopole :
        title_string{type} = 'monopole';
    elseif type == 2
        n = 6; m = 3; K = 1.3 * 10^(-2) ;	f =1;					%	Dipole :
        title_string{type} = 'dipole';
    end
    
    %	Woulff:	V_gas^n = (PxxtotW * c_air^m)/(K * rho_air* pi * R_cond * R_cond);
    V_gas{type} = ((Ea * Specvar.cel^m)/(K * f * Specvar.rho_air* pi * Specvar.R_cond^2)).^(1/n);
end
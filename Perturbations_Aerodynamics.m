function [ Acceleration_aerodynamic, Torque_aerodynamic, Force_aerodynamic, Cd_tot, Cl_tot ] = Perturbations_Aerodynamics( Spacecraft, Freestream, Surface, BODY2ECI)
%% Perturbations_Aerodynamics.m function Spacecraft, q', w', Freestream, Surface

% Purpose: Determine aerodynamic forces and torques on the spacecraft for
% given geometry, state, and freestream conditions

% Created:  Sebastian Tamrazian 2/7/2018
% Updated:  Arly Black 11/25/2019 - 

% Inputs:
%     Spacecraft: structure with fields:
%                     1. .mass : sc mass [kg]
%                     2. .inertia: spacecraft inertia [kg*m^2]
%                     3. .cm : spacecraft cm relative to geom origin [m]
%                     4. .frontal_area : spacecraft frontal area [m^2]
%                     5. .centroids : centroid locations of sc surface elements [m]
%                     6. .normals : normals of sc surface elements (unit vectors) []
%                     7. .areas : areas of sc surface elements [m^2]
%                     8. .tangents : tangents of sc surface elements (unit vectors)[]
%
%     Freestream: structure with fields:
%                     1. .density : freestream density [kg/m^3]
%                     2. .temperature : freestream temperature [K]
%                     3. .R : Freestream specific gas constant [J/(mol K)]
%                     4. .velocity : freestream velocity [m/s] 
%
%     Surface: structure with fields:
%                     1. .temperature : sc surface temperature [K]
%                     2. .reflectance_coef : sc reflectance coefficient - between 0 and 1 []
%                     3. .specular_coef : sc specular reflectance coefficient - between 0 and 1 []
%                     4. .front_emiss_coef : front emissivity coefficient - between 0 and 1 []
%                     5. .back_emiss_coef : back emissivity coefficient - between 0 and 1 []
%                     6. .front_nonLamb_coef : front non-Lambertian coefficient - between 0 and 1 []
%                     7. .back_nonLamb_coef : back non-Lambertian coefficient - between 0 and 1 []
%                     8. .transmissivity_coef = Transmissivity for use with clear sail
%                     9. .Refl                = Specify clear or reflecting surface 1=Refl, 0=Clear
%
%
%     BODY2ECI: [3x3] Rotation matrix in form [ECI] = BODY2ECI[BODY]

% Supporting Functions: 
%    FMFC - calculates free molecular flux coefficients


Freestream.beta = 1/sqrt(2*Freestream.R * Freestream.temperature);   % Beta for freestream temp
Surface.beta    = 1/sqrt(2*Freestream.R * Surface.temperature);      % Beta for wall temp (for reflected particles)
Speedratio      = Freestream.beta*norm(Freestream.velocity);         % Molecular speed ratio, using freestream temp
sigma_N         = 1;                                                 % Normal accommodation coefficient
sigma_T         = 1;                                                 % Tangential accommodation coefficient
A_ref           = Spacecraft.frontal_area;                           % Spacecraft frontal area [m^2] (drag sail area) - for use in Cd and Cl calcs


vf = BODY2ECI'*Freestream.velocity;
Vhatn = vf/norm(vf);

% Preallocate matrices for efficiency
CF_Aref = zeros(size(Spacecraft.normals))';
CM_Aref_Lref = zeros(size(Spacecraft.normals))';
Cd_A = zeros(size(Spacecraft.areas))'; 
Cl_A = zeros(size(Spacecraft.areas))';

% Iterate through all surface elements
for jj = 1:length(Spacecraft.areas)

    % Obtain element information
    Normal_jj = Spacecraft.normals(jj,:);
    Area_jj = Spacecraft.areas(jj);
    Centroid_jj = Spacecraft.centroids(jj,:);
    
    % Determine tangential component for shear calculation
    Tangential_jj = cross(Normal_jj,cross(-Vhatn,Normal_jj));
    
    % Determine angles 
    Theta = asin(dot(Vhatn,Normal_jj));     % incidence angle (angle b/w freestream velocity and tangent)
    eta = acos(dot(Vhatn,Normal_jj));       % angle offset of freestream velocity to normal = pi/2-Theta
    
    % Find pressure and shear coefficients
    [Cp,Ctau,Cd,Cl] = FMFC(sigma_N,sigma_T,Theta,Speedratio,Surface.temperature,Freestream.temperature);
   
    % Simple shadowing check - If surface is not in contact with flow, it is shadowed
    % Check to see if the surface normal is greater than 89.9 deg from the relative velocity 
    test= abs(eta);
        if test>deg2rad(89.9)
           Shadow = 0; % shadowed
        else
           Shadow = 1; % not shadowed
        end

    % Find surface cell force and moment
    CF_Aref_jj = (Cp*Normal_jj + Ctau*Tangential_jj)*Area_jj*Shadow;  % [m^2]
    CM_Aref_Lref_jj = cross(Centroid_jj, CF_Aref_jj);                 % [m^3]
      
    % Add force and moment to total
    CF_Aref(:,jj) = CF_Aref_jj;
    CM_Aref_Lref(:,jj) = CM_Aref_Lref_jj;
    
    % Find drag and lift for each surface element (Cd*A)
    Cd_A_jj = Cd*Area_jj;       % [m^2]
    Cl_A_jj = Cl*Area_jj;       % [m^2]
    
    % Add drag and lift forces to total 
    Cd_A(1,jj) = Cd_A_jj;
    Cl_A(1,jj) = Cl_A_jj;  
    
end

CF_Aref_total = sum(CF_Aref,2); % vector sum of forces [m^2]
CM_Aref_Lref_total = sum(CM_Aref_Lref,2); % vector sum of moments [m^3]
Cd_tot = sum(Cd_A,2)/A_ref; % Total spacecraft drag coefficient -> sum of Cd(i)*A(i)/Aref
Cl_tot = sum(Cl_A,2)/A_ref; % Total spacecraft lift coefficient -> sum of Cl(i)*A(i)/Aref

Freestream.q = .5 * Freestream.density * (norm(Freestream.velocity))^2;                 % Freestream dynamic pressure [N/m^2]

Force_aerodynamic = -(Freestream.q * CF_Aref_total);                                    % Force in body frame [N]
Torque_aerodynamic = -Freestream.q * CM_Aref_Lref_total;                                % Torque in body frame [Nm]
Acceleration_aerodynamic = (BODY2ECI*Force_aerodynamic)/Spacecraft.mass;                % Acceleration in ECI [m/s^2]

end

function [Cp,Ctau,Cd,Cl] = FMFC(sigma_N,sigma_T,Theta,s,T_W,T_free)
% Calculates free molecular flux coefficients

Vn = s*sin(Theta);

if isreal(Vn)
    
    Cp = (1/s^2)*(((2-sigma_N)/sqrt(pi)*Vn+sigma_N/2*sqrt(T_W/T_free))...
        *exp(-(Vn)^2)+((2-sigma_N)*((Vn)^2+1/2)+sigma_N/2*...
        sqrt(pi*T_W/T_free)*Vn)*(1+erf(Vn)));
    
    Ctau = -sigma_T*cos(Theta)/(s*sqrt(pi))*(exp(-(Vn)^2)+...
        sqrt(pi)*Vn*(1+erf(Vn)));
    
    Cd = sin(Theta)*(1+erf(Vn))*((2-sigma_N)*((sin(Theta))^2+1/(2*s^2))+...
        sigma_T*(cos(Theta))^2+sigma_N/2*sqrt(pi*T_W/T_free)*sin(Theta)/s)+...
        exp(-(Vn)^2)*((2-sigma_N)/(sqrt(pi)*s)*(sin(Theta))^2+sigma_T*(cos(Theta))^2/(sqrt(pi)*s)+...
        sigma_N/(2*s^2)*sqrt(T_W/T_free)*sin(Theta));
    
    Cl = cos(Theta)*(1+erf(Vn))*((2-sigma_N)*((sin(Theta))^2+1/(2*s^2))+...
        sigma_T*(cos(Theta))^2+sigma_N/2*sqrt(pi*T_W/T_free)*sin(Theta)/s)+...
        sin(Theta)*cos(Theta)*exp(-(Vn)^2)*((2-sigma_N)/(sqrt(pi)*s)+sigma_T/(sqrt(pi)*s)+...
        sigma_N/(2*s^2*sin(Theta))*sqrt(T_W/T_free));
    
else
    Cp=0;Ctau=0;Cd=0;Cl=0;
end

end



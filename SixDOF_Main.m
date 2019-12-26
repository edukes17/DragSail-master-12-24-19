%% SixDOF_Main.m Script

% Purpose: Simulate the six degree of freedom dynamics (orbit and attitude)
% of a spacecraft in earth orbit, taking into account environmental 
% perturbing forces and torques

% This code was written for compatibility with MATLAB 2018b
 
% Created:  Sebastian Tamrazian 10/27/2018
% Modified: Sebastian Tamrazian 3/8/2018
% Modified: Arly Black 9/24/2019
% Modified: Eileen Dukes 12/3/19

% Supporting Functions: 
%   Mesh_Data.m
%   Propagation.m
%   Convert_EA2EP.m


%% FILE PATH DEFINITIONS
% clear workspace for the simulation run
% add path for supporting folders
close all;clear;clc;
addpath(genpath('Supporting Functions'))                % Where the other functions are held, mostly conversions
addpath(genpath('Geometry Files'))                      % Where geometry files are stored

% Open file selection dialog box- choose input file
[FileName,path] = uigetfile('*.txt');
if isequal(FileName,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(path, FileName)])
end

% Identify inputs from input file
[FN,tmp,ref,spec,fem,bem,fnL,bnL,tran,Refl,mass,in,cm,aref,SMA,ecc,inc,arg,asc,TA,od,phi,theta,psi,phidot,thdot,psidot,ep,prop] = PropertiesFile(FileName);

outdir = od{1};     % Where the results will be output, you can change this to what you like, 
                    % and should change it between runs or else will overwrite previous results


%% INPUTS

%----------------Geometry----------------%

STL_Filename                = FN{1};                    % Specify filename for geometry
View_Geometry               = 1;                        % Turn on(1)/off(0) geometry viewing
View_Normals                = 1;                        % Turn on(1)/off(0) viewing surfaces normals on geometry 
View_Tangents               = 1;                        % Turn on(1)/off(0) viewing surface tangents on geometry 
View_Centroid_Labels        = 0;                        % Turn on(1)/off(0) viewing centroid labels on geometry

Scale_Geometry              = 1;                        % Note units of STL file, what to divide by to get meters (ex: STL in mm? enter 1000, in cm? enter 100) 

%---------Spacecraft properties-----------%

% "Surface"
Surface.temperature         = tmp;                      % Spacecraft surfaces temperature [K]
Surface.reflectance_coef    = ref;                      % Reflectance coef for calculating SRP effect [] non-dimensional coefficient between 0 and 1 where 1=total reflector
Surface.specular_coef       = spec;                     % Specular Reflectance  - non-dimensional coefficient between 0 and 1
Surface.front_emiss_coef    = fem;                      % Front emissivity - non-dimensional coefficient between 0 and 1
Surface.back_emiss_coef     = bem;                      % Back emissivity - non-dimensional coefficient between 0 and 1
Surface.front_nonLamb_coef  = fnL;                      % Front non-Lambertian coefficient - non-dimensional coefficient between 0 and 1
Surface.back_nonLamb_coef   = bnL;                      % Back non-Lambertian coefficient - non-dimensional coefficient between 0 and 1
Surface.transmissivity_coef = tran;                     % Transmissivity for use with clear sail
Surface.Refl                = Refl;                     % Specify clear or reflecting surface 1=Refl, 0=Clear

% "Spacecraft"
Spacecraft.mass             = mass;                     % Spacecraft mass [kg] (SRP-test1)
Spacecraft.inertia          = in;                       % Spacecraft inertia matrix [kg*m^2]
Spacecraft.cm               = cm;                       % Spacecraft center of mass, will be used to translate STL body axes to principle axes 
Spacecraft.frontal_area     = aref;                     % Dragsail frontal area [m^2]

%----------Initial orbital elements--------%

a_initial                   = SMA*1000;                 % Semimajor Axis [m]
e_initial                   = ecc;                      % Eccentricity [none]
i_initial                   = deg2rad(inc);             % Inclination [rad]
w_initial                   = deg2rad(arg);             % Argument of Perigee [rad]
o_initial                   = deg2rad(asc);             % Right Ascension of Ascending Node [rad]
t_initial                   = deg2rad(TA);              % True Anomaly [rad]

%-----------Initial orientation------------%

% Euler Angles:
phi_initial                 = deg2rad(phi);             % roll [rad] 
theta_initial               = deg2rad(theta);           % pitch [rad] 
psi_initial                 = deg2rad(psi);             % yaw [rad]

%----------Initial angular rates-----------%

phidot_initial              = deg2rad(phidot);          % roll rate [rad/s] 
thetadot_initial            = deg2rad(thdot);           % pitch rate [rad/s]
psidot_initial              = deg2rad(psidot);          % yaw rate [rad/s]

%---------Compile initial state vector-----------%

Orientation_initial         = [phi_initial,theta_initial,psi_initial];
Quaternions_initial         = Convert_EA2EP(Orientation_initial);   % convert from Euler angles to quaternions

State_initial               = [a_initial e_initial i_initial w_initial o_initial t_initial ...
                              phidot_initial thetadot_initial psidot_initial, Quaternions_initial]; 
                          
%---------------Time Settings--------------%

Time_epoch                  = ep{1};                    % Epoch, dd-mm-yyyy HH:MM:SS' [string format]

Time_propagation            = 100;                     % Duration of propagation [seconds]

Time_cell = {Time_epoch,Time_propagation};              % Contain time related parameters in cell array

         
%---------------Misc settings--------------%
Windtunnel                  = 0;                        % Specify 1 for Windtunnel/3DOF mode, 0 for Orbital/6DOF mode

Tolerance                   = 1e-6;                     % Specify ODE solver tolerance

global Deorb_Alt;                                       % Global variable- defined for all associated functions
Deorb_Alt = 110*1000;                                   % Deorbit altitude [m]


%% OBTAIN GEOMETRY INFORMATION
[Spacecraft.centroids, Spacecraft.normals, Spacecraft.areas, ~, Spacecraft.tangents] = MeshData(STL_Filename, View_Geometry, Scale_Geometry, Spacecraft, View_Normals, View_Tangents, View_Centroid_Labels, outdir);

%% PROPAGATE SPACECRAFT MOTION

[ Time, State ] = Propagation( Windtunnel, State_initial, Orientation_initial, Time_cell, Tolerance, Deorb_Alt, Spacecraft, Surface, outdir);


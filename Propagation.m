function [ Time_history, State_history ] = Propagation( Windtunnel, State_initial, Orientation_initial, Time_cell, Tolerance, Deorb_Alt, Spacecraft, Surface, outdir, OS)
%% Propagation.m function

%Purpose: Propagate the attitude and trajectory of the spacecraft in orbit.
 
%Created:  Sebastian Tamrazian 12/23/2018
%Modified: Sebastian Tamrazian 2/12/2019 
%Modified: Arly Black 9/24/2019
 
%Acknowledgement: Jay Iuliano for his excellent help and knowledge of orbit
%propagation methods

%Inputs:
%     Windtunnel: 1 or 0 to specify if the mode is active/inactive
%
%     State_initial: Initial user defined state of the spacecraft [1x12]
%                  Initial set of classical orbital elements in the order:
%                     1. Semimajor Axis (m)
%                     2. Eccentricity (-)
%                     3. Inclination (rad)
%                     4. Argument of Perigee (rad)
%                     5. Right Ascension of the Ascending Node (rad)
%                     6. True Anomaly (rad)
%                  Initial Euler angles are in the order:
%                     1. Phi (rad)
%                     2. Theta (rad)
%                     3. Psi (rad)
%                     4. Phidot (rad/s)
%                     5. Thetadot (rad/s)
%                     6. Psidot (rad/s)
%
%     Orientation_initial: [3x1] vector specifying initial euler angle
%     state
%
%     Time_cell: User defined time options for propagation {1x3}
%                     1. Initial date: 'dd-mm-yyyy HH:MM:SS' (string)
%                     2. Propagation time (seconds)
%
%     Tolerance: Specifies solver numerical tolerance 
%
%     Deorb_Alt: User specified deorbit altitude (default 110 km) [m]
%     
%     Spacecraft: structure with fields:
%                     1. .mass : sc mass [kg]
%                     2. .inertia: spacecraft inertia [kg*m^2]
%                     3. .cm : spacecraft cm relative to geom origin [m]
%                     4. .centroids : centroid locations of sc surface elements [m]
%                     5. .normals : normals of sc surface elements (unit vectors) []
%                     6. .areas : areas of sc surface elements [m^2]
%                     7. .tangents : tangents of sc surface elements (unit vectors)[]
%
%     Surface: structure with fields:
%                     1. .temperature : sc surface temperature [K]
%                     2. .reflectance_coef : sc reflectance coefficient - between 0 and 1 []
%                     3. .specular_coef : sc specular reflectance coefficient - between 0 and 1 []
%                     4. .front_emiss_coef : front emissivity coefficient - between 0 and 1 []
%                     5. .back_emiss_coef : back emissivity coefficient - between 0 and 1 []
%                     6. .front_nonLamb_coef : front non-Lambertian coefficient - between 0 and 1 []
%                     7. .back_nonLamb_coef : back non-Lambertian coefficient - between 0 and 1 []
%
%     outdir: specifies output file path for simulation results
%
%     OS: specifies operating system type: PC or Mac (impacts file directory storage syntax) 


% Supporting Functions: 
%    Perturbations_Aerodynamics.m
%    Perturbations_Gravitational.m
%    Perturbations_SolarRadiationPressure.m
%    AtmosphericProperties.m
%    eclipseCheck.m
%    SunPosition.m
%    Convert_BODY2ECI.m
%    Convert_COES2ECI.m
%    Convert_EA2EP.m
%    PlotResults.m
                                  
%% Initialize Variables

% Earth Properties
EARTH.SGP                   = 3.9860044189e14;          % Std. grav. param. [m^3/s^2]
EARTH.MASS                  = 5.97e24;                  % Mass [kg]
EARTH.EQRADIUS              = 6378.1363e3;              % Equatorial radius [m]
EARTH.PORADIUS              = 6051.8e3;                 % Polar radius [m]
EARTH.J2                    = 1082.63e-6;               % J2 oblateness parameter [-]


% Initial orbital elements (classical)
a_initial                   = State_initial(1);         % Semi-major axis [m] 
e_initial                   = State_initial(2);         % Eccentricity []
i_initial                   = State_initial(3);         % Inclination [rad]
w_initial                   = State_initial(4);         % Arguement of periapsis [rad]
o_initial                   = State_initial(5);         % Raan [rad]
t_initial                   = State_initial(6);         % True anomaly [rad]

% Initial attitude state
Angularrates_initial        = State_initial(7:9);       % [rad/s]
Quaternions_initial         = State_initial(10:13);     % [q1 q2 q3 q4]

% Initial radius and velocity vectors from classical orbital elements
AngularMomentum_initial = sqrt(a_initial*(1-e_initial^2)*EARTH.SGP);    % angular momentum [m2/s]
[Position_initial,Velocity_initial] = Convert_COES2ECI(AngularMomentum_initial,t_initial,w_initial,o_initial,i_initial,e_initial,EARTH.SGP); 


% Initial state vectors: SixDOF_initial is the intial state of the
% integrator
% Replace complex transpose ' with real transpose .'
ThreeDOF_initial            = [Quaternions_initial,Angularrates_initial]; % Windtunnel mode
SixDOF_initial              = [Position_initial.', Velocity_initial.', Quaternions_initial, Angularrates_initial]; % Orbital mode

% Time Options
Time_epoch                  = Time_cell{1};    % 'mm-dd-yyyy HH:MM:SS'
Time_propagation            = Time_cell{2};    % [seconds]

%% Propagate 3DOF or 6DOF depending on if Windunnel mode is active

if ~Windtunnel 
    %Runs if Orbital mode (6DOF) is enabled
    
    % Set options for specified tolerance and check for deorbit
    options  = odeset('Abstol', Tolerance, 'Reltol', Tolerance, 'Events', @Deorbit);
    
    % Find position of sun
    [R_sun] = SunPosition(Time_epoch);
 
    %define timeseries to save parameters    
   
    global tsAll
    
    tsAll = timeseries('AllStuff');
    
    profile on 
   
    %Propagate
    [Time_history,State_history] = ode113(@(Time, State) Orbitalcase(Time, State, Spacecraft, Surface, EARTH, Time_epoch, R_sun),...
                                            [0 Time_propagation], SixDOF_initial, options);
    profile viewer

    State_history = real(State_history);
    
    %% Output results
    
    [T1,T2,T3,T4,T5,T6,T7,T8]=UnpackSavedData(tsAll,Time_propagation);
   
        
    % Check if user OS is PC or Mac, then output results to indicated folder
    if ispc 
        % For PC
      writetable(T1,[pwd '\' outdir '\AtmosAngles.csv']);
      writetable(T2,[pwd '\' outdir '\Euler.csv']);
      writetable(T3,[pwd '\' outdir '\Aero.csv']);
      writetable(T4,[pwd '\' outdir '\Coeffs.csv']);
      writetable(T5,[pwd '\' outdir '\Solar.csv']);
      writetable(T6,[pwd '\' outdir '\Grav.csv']);
%      writetable(T7,[pwd '\' outdir '\All.csv']);
      writetable(T8,[pwd '\' outdir '\TotalForces.csv']);
    elseif ismac
        % For Mac
        writetable(T1,[pwd '/' outdir '/AtmosAngles.csv']);
        writetable(T2,[pwd '/' outdir '/Euler.csv']);
        writetable(T3,[pwd '/' outdir '/Aero.csv']);
        writetable(T4,[pwd '/' outdir '/Coeffs.csv']);    
        writetable(T5,[pwd '/' outdir '/Solar.csv']);
        writetable(T6,[pwd '/' outdir '/Grav.csv']);
        writetable(T7,[pwd '/' outdir '/All.csv']);
        writetable(T8,[pwd '/' outdir '/TotalForces.csv']);
    end
    
    % Plot results
    PlotResults(Time_history, State_history,R_sun,Orientation_initial,EARTH,T1,T2,T3,T4,T5,T6,T8,outdir);

end
end


%% Orbital Formulation ODE function

 function [Statederivative] = Orbitalcase(Time, State, Spacecraft, Surface, EARTH, Time_epoch, R_sun)
          
   global tsAll
   
 % Translational state
    r = State(1:3); %position [m]  
    v = State(4:6); %velocity [m/s]  
 
    % Rotational state
    q = real(State(7:10));                          % quaternions, q(4) is the scalar component
    w = State(11:13);                               % angular rates [rad/s]
    Euler_angles = Convert_EP2EA(q(1:4).');
    [BODY2ECI] = Convert_BODY2ECI(Euler_angles);    % body to ECI coordinate conversion
    Euler = rad2deg(Euler_angles);
        
    % Obtain freestream conditions and earth shadow for current state
    [Freestream] = AtmosphericProperties(EARTH, Time_epoch, r, v);
    vr_body = BODY2ECI.'*Freestream.velocity;
    vrhat = vr_body/norm(vr_body);
    %Based on Arun's defn, alpha is the X/Y and beta is in the X/Z
    alpha = atan2(vrhat(2),vrhat(1));             % calculate angle of attack [rad]
    ang_of_attack = rad2deg(alpha.');
    beta = atan2(vrhat(3),vrhat(1));              % calculate side slip angle [rad]
    sideslip = rad2deg(beta);
    vr_out = transpose(vrhat);
    aoa = acos(dot(vrhat,[1,0,0]));              % total angle of attack [rad] (angle bw vr and x axis)
    aoa_tot=rad2deg(aoa);
        
    [Shadow] = eclipseCheck(R_sun,r,EARTH);
    
    % Find  perturbations
    [Accel_aerodynamic,Torque_aerodynamic,Force_aerodynamic, Cd_tot, Cl_tot] = Perturbations_Aerodynamics(Spacecraft, Freestream, Surface, BODY2ECI);
    
    [Accel_aspherical, Torque_gravity,~] = Perturbations_Gravitational(r, EARTH, Spacecraft, BODY2ECI);
    
    [Accel_srp, Torque_srp,Force_srp,R_sun_hat] = Perturbations_SolarRadiationPressure(Spacecraft, r, R_sun, Surface, Shadow, BODY2ECI);
        
    % Sum Perturbations
    Accelerations = Accel_aerodynamic + Accel_aspherical + Accel_srp;
    Torques = Torque_aerodynamic + Torque_gravity + Torque_srp; %Torques should be a vector, was a 3x3

  % Translation :
    a = -EARTH.SGP/norm(r)^3*r + Accelerations;
       
  % Rotation :
    % Dynamic EOMs:
    % Torque = I*wdot + (w x H)
    H = Spacecraft.inertia*w;           % angular momentum = inertia tensor * angular velocity [kg*m2/s]
    BKE = cross(w,H);                   % w x H [kg*m2/s2 = Nm]
    I_wdot = Torques-BKE;               % I*wdot = Torque - (w x H) [Nm]
    wd = I_wdot./Spacecraft.inertia;    % wdot matrix
    wdot = [wd(1,1);wd(2,2);wd(3,3)];   % angular rate derivatives [rad/s2]

    % Kinematic EOMs (quaternion formulation):
    qdot(1,1) = (w(1)*q(4)-w(2)*q(3)+w(3)*q(2))/2;
    qdot(2,1) = (w(1)*q(3)+w(2)*q(4)-w(3)*q(1))/2;
    qdot(3,1) =(-w(1)*q(2)+w(2)*q(1)+w(3)*q(4))/2;
    qdot(4,1) =-(w(1)*q(1)+w(2)*q(2)+w(3)*q(3))/2;

    % 6DOF state derivative
    Statederivative = [v;a;qdot;wdot];

    %   Save all of the variables 
    Euler_out = transpose(Euler);
%    w_out = transpose(w);
    SaveData = [ang_of_attack;sideslip;aoa_tot;vrhat;Torque_aerodynamic;Force_aerodynamic;Cd_tot;Cl_tot;Torque_gravity;Accel_aspherical;Torque_srp;Force_srp;R_sun_hat;Torques;Accelerations;Euler_out;w];
    tsAll = addsample(tsAll,'Time',Time,'OverwriteFlag',true,'Data',SaveData.');
 
    % Print Status
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b') % delete previous time
    fprintf('%13.2f',Time) % print current time

    
 end

function [value,isterminal,direction] = Deorbit(~, State, ~, ~, ~, ~, ~, ~)    

    r = State(1:3); % position [m]
    global Deorb_Alt
    
    %determine altitude
    Altitude = norm(r)-6378100;
    if Altitude <= Deorb_Alt
        check = 1;
    else
        check = 0;
    end
    
    %if Altitude is less than "Deorb_Alt" m, propagation stops
    value = 1 - check;
    isterminal = 1;
    direction = 0;
end
 
function [T1,T2,T3,T4,T5,T6,T7,T8] = UnpackSavedData(tsAll,Time_propagation)

% SaveData =
% [ang_of_attack;sideslip;aoa_tot;vrhat;Torque_aerodynamic;Force_aerodynamic;Cd_tot;...
% Cl_tot;Torque_gravity;Accel_aspherical;Torque_srp;Force_srp;R_sun_hat;Torques;Accelerations;Euler,w]
% Convert timeseries to a table  
    T7 = table(tsAll.time,tsAll.data,'VariableNames',{'Time','SaveData'});

% Unpack SaveData to variables  
    Saved_Data = T7.('SaveData');
    Alpha_out = Saved_Data(:,1);
    Beta_out = Saved_Data(:,2);
    AOA_out = Saved_Data(:,3);
    vr_out_x = Saved_Data(:,4);
    vr_out_y = Saved_Data(:,5);
    vr_out_z = Saved_Data(:,6);
    Taero_out_x = Saved_Data(:,7);
    Taero_out_y = Saved_Data(:,8);
    Taero_out_z = Saved_Data(:,9);
    Faero_out_x = Saved_Data(:,10);
    Faero_out_y = Saved_Data(:,11);
    Faero_out_z = Saved_Data(:,12);
    Cd_out = Saved_Data(:,13);
    Cl_out = Saved_Data(:,14);
    Tgrav_out_x = Saved_Data(:,15);
    Tgrav_out_y = Saved_Data(:,16);
    Tgrav_out_z = Saved_Data(:,17);
    ACgrav_out_x = Saved_Data(:,18);
    ACgrav_out_y = Saved_Data(:,19);
    ACgrav_out_z = Saved_Data(:,20);
    Tsol_out_x = Saved_Data(:,21);
    Tsol_out_y = Saved_Data(:,22);
    Tsol_out_z = Saved_Data(:,23);
    Fsol_out_x = Saved_Data(:,24);
    Fsol_out_y = Saved_Data(:,25);
    Fsol_out_z = Saved_Data(:,26);
    Sbody_out_x = Saved_Data(:,27);
    Sbody_out_y = Saved_Data(:,28);
    Sbody_out_z = Saved_Data(:,29);
    Tork_out_x = Saved_Data(:,30);
    Tork_out_y = Saved_Data(:,31);
    Tork_out_z = Saved_Data(:,32);
    Accel_out_x = Saved_Data(:,33);
    Accel_out_y = Saved_Data(:,34);
    Accel_out_z = Saved_Data(:,35);
    Euler_out_z = Saved_Data(:,36);
    Euler_out_y = Saved_Data(:,37);
    Euler_out_x = Saved_Data(:,38);
    w_out_x = Saved_Data(:,39);
    w_out_y = Saved_Data(:,40);
    w_out_z = Saved_Data(:,41);
    
 % Make output tables   
    T1 = table(T7.('Time'),Alpha_out,Beta_out,AOA_out,'VariableNames',{'Time','Alpha','Beta','AOA'});
    T2 = table(T7.('Time'),Euler_out_z,Euler_out_y,Euler_out_x,w_out_x,w_out_y,w_out_z,'VariableNames',...
        {'Time','Roll','Pitch','Yaw','BodyRateX','BodyRateY','BodyRateZ'});
    T3 = table(T7.('Time'),Taero_out_x,Taero_out_y,Taero_out_z,Faero_out_x,Faero_out_y,Faero_out_z,...
        'VariableNames',{'Time','AeroTx','AeroTy','AeroTz','AeroFx','AeroFy','AeroFz'});
    T4 = table(T7.('Time'),Cd_out,'VariableNames',{'Time','Cd'});
% Set the time interval for calculating the average Cd    
        if Time_propagation > 1000      
          k = 1000;
        else
           k=10; 
        end   
    Cd_ave = movmean(T4.('Cd'),k,'SamplePoints',T4.('Time'));
    T4 = table(T7.('Time'),Cd_out,Cd_ave,Cl_out,'VariableNames',{'Time','Cd','Cd_ave','Cl'});
    T5 = table(T7.('Time'),Tsol_out_x,Tsol_out_y,Tsol_out_z,Fsol_out_x,Fsol_out_y,Fsol_out_z,Sbody_out_x,...
        Sbody_out_y,Sbody_out_z,'VariableNames',{'Time','SolarTx','SolarTy','SolarTz','SolarFx','SolarFy',...
        'SolarFz','SunBodyX','SunBodyY','SunBodyZ'});
    T6 = table(T7.('Time'),Tgrav_out_x,Tgrav_out_y,Tgrav_out_z,ACgrav_out_x,ACgrav_out_y,...
        ACgrav_out_z,'VariableNames',{'Time','GravTx','GravTy','GravTz','GravAx','GravAy','GravAz'});
    T8 = table(T7.('Time'),Tork_out_x,Tork_out_y,Tork_out_z,Accel_out_x,Accel_out_y,Accel_out_z,...
        'VariableNames',{'Time','Torque_X', 'Torque_Y','Torque_Z','Accel_X','Accel_Z','Accel_Y'});
   
    
end    

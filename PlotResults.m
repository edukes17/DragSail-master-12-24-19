function PlotResults(Time_history, State_history,R_sun,~,EARTH,T2,T3,T4,outdir)

if ~exist(outdir,'dir') %make output directory if it doesn't exist already
    mkdir(outdir)
end
 global ts1
 global ts2
 global ts3
 global ts4
 global ts5
 global ts6
 global ts7
 global ts8
 global tsc9
 global ts9
 global tsf1
 global tsf2
 global tsf3
 global tsf4
 global tst1


%---------Cd and Cl History------------------------%
figure(2)
coeff_drag = subplot(2,1,1);

hold on
plot(ts7);
plot(T4.('Time'),T4.('Cd_ave'),'m','Linewidth',2);
hold off
legend('Cd','Cd ave');
ylabel('Coeff Drag');xlabel('Time [sec]');
title(coeff_drag, 'Cd history');

coeff_lift = subplot(2,1,2);
plot(ts8);
legend('Cl');
ylabel('Coeff Lift');xlabel('Time [sec]');
title(coeff_lift,'Cl history');

if ispc
    savefig([pwd '\' outdir '\' 'Coeffs'])
elseif ismac
    saveas(gcf,[pwd '/' outdir '/' 'Coeffs.jpg'])
end
 
 
%---------------------Euler Angle History---------------------%
figure(3)       
Rollplot = subplot(3,1,1);
Angles = T2.('Euler');
Roll = Angles(:,1);
plot(T2.('Time'),Roll,'m','LineWidth',1);
legend('\phi (roll)');
ylabel('Angle [deg]');xlabel('Time [sec]');
title(Rollplot, 'Roll Angle');

pitchplot = subplot(3,1,2);
Pitch = Angles(:,2);
plot(T2.('Time'),Pitch,'g','LineWidth',1);
legend('\theta (pitch)');
ylabel('Angle [deg]');xlabel('Time [sec]');
title(pitchplot, 'Pitch Angle');

yawplot = subplot(3,1,3);
Yaw = Angles(:,3);
plot(T2.('Time'),Yaw,'LineWidth',1);
legend('\psi (yaw)');
ylabel('Angle [deg]');xlabel('Time [sec]');
title(yawplot, 'Yaw Angle');

if ispc
    savefig([pwd '\' outdir '\' 'eulerplots'])
elseif ismac
    saveas(gcf,[pwd '/' outdir '/' 'eulerplots.jpg'])   %changed to saveas

end


%---------------------Altitude History---------------------%
figure(3)
Altitude_history = zeros(1,length(State_history));
for nn = 1:length(State_history)
    Altitude_history(nn) = norm(State_history(nn,1:3))/1000-EARTH.EQRADIUS/1000;
    %         Altitude_history(nn) = norm(State_history(nn,1:3))/1000;
end
plot(Time_history,Altitude_history);
xlabel('Time [sec]');ylabel('Altitude [km]');
title('Altitude history');

if ispc
    savefig([pwd '\' outdir '\' 'altplot'])
elseif ismac
    saveas(gcf, [pwd '/' outdir '/' 'altplot.jpg'])
end


%---------------------Aero torques and forces---------------------%

figure(5)
torqueplotaeroX = subplot(5,1,1);
AeroT = T3.('AeroT');
AeroTx = AeroT(:,1);
plot(T3.('Time'),AeroTx);
%legend('AeroTx');

ylabel('Torque [Nm]');xlabel('Time [sec]');
title(torqueplotaeroX, 'X Aero torque history');

torqueplotaeroY = subplot(5,1,2);
AeroTy = AeroT(:,2);
%legend('AeroTy');
ylabel('Torque [Nm]');xlabel('Time [sec]');
title(torqueplotaeroY, 'Y Aero torque history');

torqueplotaeroZ = subplot(5,1,3);
AeroTz = AeroT(:,3);
plot(T3.('Time'),AeroTz);
%legend('AeroTz');
ylabel('Torque [Nm]');xlabel('Time [sec]');
title(torqueplotaeroZ, 'Z Aero torque history');

%forceplotaero = subplot(5,1,1,2);
%plot(ts2);
%legend('Fx','Fy','Fz');
%ylabel('Force [N]');xlabel('Time [sec]');
%title(forceplotaero,'Aero force history');

if ispc
    savefig([pwd '\' outdir '\' 'aeroplots'])
elseif ismac
    saveas(gcf,[pwd '/' outdir '/' 'aeroplots.jpg'])
end


%---------------------Gravity torques and forces---------------------%

figure(6)
plot(ts3);

legend('Tx','Ty','Tz');
ylabel('Torque [Nm]');xlabel('Time [sec]');
title('Gravity gradient torque history');

%forceplotgrav = subplot(2,1,2);
%hold on
%plot(Time_history,Forces_grav(1,:));
%plot(Time_history,Forces_grav(2,:));
%plot(Time_history,Forces_grav(3,:));
%hold off
%legend('Fx','Fy','Fz');
%ylabel('Force [N]');xlabel('Time [sec]');
%title(forceplotgrav,'Gravity force history');

if ispc
    savefig([pwd '\' outdir '\' 'gravplots'])

elseif ismac
    saveas(gcf,[pwd '/' outdir '/' 'gravplots.jpg'])

end


%-------------------SRP torques and forces-----------------------%
figure(6)
torqueplotsrp = subplot(2,1,1);
plot(ts4);
legend('Tx','Ty','Tz');
ylabel('Torque [Nm]');xlabel('Time [sec]');
title(torqueplotsrp, 'SRP torque history');

forceplotsrp = subplot(2,1,2);
plot(ts5);
legend('Fx','Fy','Fz');
ylabel('Force [N]');xlabel('Time [sec]');
title(forceplotsrp,'SRP force history');

if ispc
    savefig([pwd '\' outdir '\' 'srpplots'])
elseif ismac
    saveas(gcf,[pwd '/' outdir '/' 'srpplots.jpg'])

end


%--------Angle of attack, Side Slip Angle, Total Angle of Attack---------%
figure(7)
alphaplot = subplot(3,1,1);
plot(tsf1);
xlabel('Time [sec]');ylabel(['\alpha [deg]']);
% ylabel('$\alpha$ [deg]','Interpreter','latex');
title(alphaplot, 'Angle of attack history');

betaplot = subplot(3,1,2);
plot(tsf2);
xlabel('Time [sec]');ylabel(['\beta [deg]']);
% ylabel('$\beta$ [deg]','Interpreter','latex');
title(betaplot, 'Side slip angle history');

aoaplot = subplot(3,1,3);
plot(tsf3);
xlabel('Time [sec]');ylabel(['\alpha_T [deg]']);
title(aoaplot, 'Total Angle of Attack history');


if ispc
    savefig([pwd '\' outdir '\' 'FlowAngleplots'])
elseif ismac
    saveas(gcf,[pwd '/' outdir '/' 'FlowAngleplots.jpg'])

end


%-------------------Orbit and Sun Direction----------------------%
figure(8)
hold on
eqrad = EARTH.EQRADIUS/1000;
porad = EARTH.PORADIUS/1000;
rsun = R_sun/1000;
rx = State_history(:,1)/1000;
ry = State_history(:,2)/1000;
rz = State_history(:,3)/1000;
[Earthplotx,Earthploty,Earthplotz] = ellipsoid(0,0,0,eqrad,eqrad,porad);
surf(Earthplotx,Earthploty,Earthplotz,'FaceColor', [0 .5 .5],'EdgeColor',[0 .7 0]);
alpha 0.7  % transparency of surface plot (0 = fully transparent)
plot3(rx,ry,rz,'Linewidth',.8,'Color','k');
xlabel('x[km]');ylabel('y[km]');zlabel('z[km]');
Sun_vec = rsun/norm(rsun);
quiver3(0,0,0,Sun_vec(1)*9500,Sun_vec(2)*9500,Sun_vec(3)*9500,'LineWidth',2,'MaxHeadSize',0.7,'Color','m');
legend('Earth','Orbit path','Sun direction');
axis equal
axis auto
hold off
title('Orbit');

if ispc
    savefig([pwd '\' outdir '\' 'orbitplot'])

elseif ismac
    saveas(gcf,[pwd '/' outdir '/' 'orbitplot.jpg'])
end

ylabel('Sun unit vector (body frame)');xlabel('Time [sec]')
title(['Sun Unit Vector in Body Frame'])

if ispc
    savefig([pwd '\' outdir '\' 'sunvect'])

elseif ismac
    saveas(gcf,[pwd '/' outdir '/' 'sunvect.jpg'])
end

%-----------Total Torques------------------------%
figure(11)
plot(tst1,'LineWidth',2); %Total
xlabel('Time [sec]');ylabel('Torque');
legend('TTx','TTy','TTz');
title('Total Torque');

if ispc
    savefig([pwd '\' outdir '\' 'Torques'])
elseif ismac
    saveas(gcf,[pwd '/' outdir '/' 'Torques.jpg'])
end

%---------------Total Angle of Attack--------------------%
figure(10)
plot(Time_history,rad2deg(aoa_tot));
xlabel('Time [sec]');ylabel(['\alpha_T [deg]']);
title('Total Angle of Attack history');

if strcmp(OS,'PC') == 1
    savefig([pwd '\' outdir '\' 'aoaplot'])
elseif strcmp(OS,'Mac') == 1
    savefig([pwd '/' outdir '/' 'aoaplot'])
end

%-----------Cd vs alpha vs beta surface plot-------------%
figure(11)
plot3(rad2deg(aoa),rad2deg(beta),Cd,'x');
xlabel(['\alpha [deg]']);ylabel(['\beta [deg]']);zlabel(['C_d']);
xlim([0 180]);ylim([0 180]);


% figure(12)
% Cdalphaplot = subplot(2,1,1);
% hold on
% plot(rad2deg(aoa),Cd);
% xlabel(['\alpha [deg]']);ylabel(['C_d']);
% % ylabel('$\alpha$ [deg]','Interpreter','latex');
% title(Cdalphaplot, 'Cd vs Angle of Attack');
% 
% Cdbetaplot = subplot(2,1,2);
% hold on
% plot(rad2deg(beta),Cd);
% xlabel(['\beta [deg]']);ylabel(['C_d']);
% ylabel('$\beta$ [deg]','Interpreter','latex');
% title(Cdbetaplot, 'Cd vs Side Slip Angle');

end
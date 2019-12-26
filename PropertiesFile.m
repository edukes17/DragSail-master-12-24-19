function [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,bb,cc] = PropertiesFile(FN)    

% Enter filename
FileName = sprintf(FN); 

% Open file
fid = fopen(FileName);
matData = textscan(fid,'%s%s%s%s%s%s%s', 'Delimiter',{',','=','%'}, 'CollectOutput',true);

% Identify rows in which desired constants are located
idx1 = find(contains(matData{1}, 'STL Filename'),1);     % Search for text containing "STL Filename" and identify the row number
idx2 = find(contains(matData{1}, 'temperature'),1);         %b
idx3 = find(contains(matData{1}, 'reflectance_coef'),1);    %c
idx4 = find(contains(matData{1}, 'specular_coef'),1);       %d
idx5 = find(contains(matData{1}, 'front_emiss_coef'),1);    %e
idx6 = find(contains(matData{1}, 'back_emiss_coef'),1);     %f
idx7 = find(contains(matData{1}, 'front_nonLamb_coef'),1);  %g
idx8 = find(contains(matData{1}, 'back_nonLamb_coef'),1);   %h
idx9 = find(contains(matData{1}, 'transmissivity_coef'),1); %i
idx10 = find(contains(matData{1}, 'Refl'),1);               %j
idx11 = find(contains(matData{1}, 'mass'),1);               %k
idx12 = find(contains(matData{1}, 'inertia'),1);            %l
idx13 = find(contains(matData{1}, 'cm'),1);                 %m
idx14 = find(contains(matData{1}, 'frontal_area'),1);       %n
idx15 = find(contains(matData{1}, 'Semimajor'),1);          %o
idx16 = find(contains(matData{1}, 'Eccentricity'),1);       %p
idx17 = find(contains(matData{1}, 'Inclination'),1);        %q
idx18 = find(contains(matData{1}, 'Argument'),1);           %r
idx19 = find(contains(matData{1}, 'Ascension'),1);          %s
idx20 = find(contains(matData{1}, 'True Anomaly'),1);       %t
idx21 = find(contains(matData{1}, 'Output Directory'),1);   %u
idx22 = find(contains(matData{1}, 'Roll angle'),1);         %v
idx23 = find(contains(matData{1}, 'Pitch angle'),1);        %w
idx24 = find(contains(matData{1}, 'Yaw angle'),1);          %x
idx25 = find(contains(matData{1}, 'Roll rate'),1);          %y
idx26 = find(contains(matData{1}, 'Pitch rate'),1);         %z
idx27 = find(contains(matData{1}, 'Yaw rate'),1);           %aa
idx28 = find(contains(matData{1}, 'Epoch'),1);              %bb
idx29 = find(contains(matData{1}, 'Propagation'),1);        %cc

fclose(fid);

info = matData{1,1};

% Identify which row and column the relevant value is in. Convert from
% string to double
a = (info(idx1,2));             %filename
b = str2double(info(idx2,2));   %temperature
c = str2double(info(idx3,2));   %refl_c
d = str2double(info(idx4,2));
e = str2double(info(idx5,2));
f = str2double(info(idx6,2));
g = str2double(info(idx7,2));
h = str2double(info(idx8,2));
i = str2double(info(idx9,2));
j = str2double(info(idx10,2));  %refl,yes/no
k = str2double(info(idx11,2));  %mass
l = [str2double(info(idx12,2)),str2double(info(idx12,3)),str2double(info(idx12,4));str2double(info(idx12+1,2)),str2double(info(idx12+1,3)),str2double(info(idx12+1,4));str2double(info(idx12+2,2)),str2double(info(idx12+2,3)),str2double(info(idx12+2,4))];
m = [str2double(info(idx13,2)),str2double(info(idx13,3)),str2double(info(idx13,4))];
n = str2double(info(idx14,2));  %frontal area
o = str2double(info(idx15,2));  %semimajor
p = str2double(info(idx16,2));  %ecc
q = str2double(info(idx17,2));  %inc
r = str2double(info(idx18,2));  %arg
s = str2double(info(idx19,2));  %asc
t = str2double(info(idx20,2));  %true
u = (info(idx21,2));            %outputdir
v = str2double(info(idx22,2));  %roll
w = str2double(info(idx23,2));  %pitch
x = str2double(info(idx24,2));  %yaw
y = str2double(info(idx25,2));  %roll rate
z = str2double(info(idx26,2));  %pitch rate
aa = str2double(info(idx27,2)); %yaw rate
bb =(info(idx28,2));            %epoch
cc = str2double(info(idx29,2)); %prop

end
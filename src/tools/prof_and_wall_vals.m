% extract_perchan_profile.m
%
% Main program to extract the turbulent profile:
function prof_and_wall_vals(geom,model,wallfile)

% Extract data
[XY,scalars,turb1,turb2,turb3,turb4,muLam,dudy,muT,uu,vv,ww,uv] ...
    = importSSTprofile(model);

[rho,vel,press] = scalar_decomp(scalars);

% Summarize Data about current flow 
srho   = rho(1,end);
spress = press(1,end);
svel   = vel(1,end,1);
smulam = muLam(1,end);

freestream = [svel, spress, srho, smulam];

fprintf('Free stream density = %.4f\n',srho);
fprintf('Free stream pressure = %.2f\n',spress);
fprintf('Free stream velocity = %.3f\n',svel);
fprintf('Free stream laminar viscosity = %.3e\n\n',smulam);

% Calculate layer thickness
[delta delta_idx] = blthickness(XY,vel);

% Calculate Re_theta and Re_disp
[dispth,theta,Re_disp,Re_theta] = calc_Retheta(XY,vel,rho,muLam,freestream,delta_idx);

% Calculate wall profiles and determine location for extraction
if strcmp(geom,'channel')
    [Re_tau] = gen_wall_values(wallfile,freestream,Re_theta,ones(1,length(delta)));
    [imatch fak] = extract_point(XY(:,1,1));
elseif strcmp(geom,'plate')
    [Re_tau] = gen_wall_values(wallfile,freestream,Re_theta,delta);
    [imatch fak] = extract_point(Re_theta);
else
    fprintf('Unknown geometry specified\n');
end

% Generate cross-sectional profile
gen_profile(XY,vel,rho,press,turb1,turb2,turb3,turb4,muLam,muT,...
    uu,vv,ww,uv,imatch,fak,model);

% Plot the variation of various quantities through the provided data,
% illustrating where the turbulent profile has been extracted
close all;
val_delta    = delta(imatch-1)    + fak*(delta(imatch)    - delta(imatch-1)   );
val_dispth   = dispth(imatch-1)   + fak*(dispth(imatch)   - dispth(imatch-1)  );
val_theta    = theta(imatch-1)    + fak*(theta(imatch)    - theta(imatch-1)   );
val_Re_disp  = Re_disp(imatch-1)  + fak*(Re_disp(imatch)  - Re_disp(imatch-1) );
val_Re_theta = Re_theta(imatch-1) + fak*(Re_theta(imatch) - Re_theta(imatch-1));
val_Re_tau   = Re_tau(imatch-1)   + fak*(Re_tau(imatch)   - Re_tau(imatch-1)  );

figure(1)
subplot(2,2,1);
plot(XY(:,1,1),delta(:,1,1),'-o',[XY(imatch,1,1) XY(imatch,1,1)],[val_delta*.8 val_delta*1.2],'--r');
title('BL height (99%)');xlabel('X (m)');

subplot(2,2,2);
plot(XY(:,1,1),Re_disp(:),'-o',[XY(imatch,1,1) XY(imatch,1,1)],[val_Re_disp*.95 val_Re_disp*1.05],'--r');
title('Re_{disp}');xlabel('X (m)');

subplot(2,2,3);
plot(XY(:,1,1),Re_tau(:),'-o',[XY(imatch,1,1) XY(imatch,1,1)],[val_Re_tau*.8 val_Re_tau*1.2],'--r');
title('Re_{tau}');xlabel('X (m)');

subplot(2,2,4);
plot(XY(:,1,1),Re_theta(:),'-o',[XY(imatch,1,1) XY(imatch,1,1)],[val_Re_theta*.8 val_Re_theta*1.2],'--r');
title('Re_{theta}');xlabel('X (m)');

fprintf('\nBL height  : %.5f\n',val_delta);
fprintf('Disp thick : %.5f\n',val_dispth);
fprintf('Mom thick  : %.5f\n\n',val_theta);
fprintf('Re_disp    : %.3f\n',val_Re_disp);
fprintf('Re_theta   : %.3f\n',val_Re_theta);
fprintf('Re_tau     : %.3f\n',val_Re_tau);
end

% ------------------------ FUNCTION FILES ---------------------------------

% Function to import the results from the flat-plate simulation
function [XY,scalars,turb1,turb2,turb3,turb4,muLam,dudy,muT,uu,vv,ww,uv]...
    = importSSTprofile(model)
% This reads in unstructured output profile, then sorts by X and Y location
% to get into structured (i,j) format

file = 'flowfield.txt';
Data = importdata(file);

% Spatial locations
XY = Data(:,1:2);

% Compute imax (dimension in y direction)
sortedXY = sortrows(XY,[1 2]);
imax = 1;
val = sortedXY(1,2);
while sortedXY(imax+1,2) > val
    imax = imax+1;
end

% Compute jmax (dimension in x direction)
jmax = size(Data,1) / imax;

% Sort in X
Data = sortrows(Data,[1]);

% Spatial locations
XY = zeros(jmax,imax,2);

% Regular Scalars - rho,Xvel,Yvel,press
scalars = zeros(jmax,imax,4);

% Turbulence model scalars
turb1 = zeros(jmax,imax);
turb2 = zeros(jmax,imax);
turb3 = zeros(jmax,imax);
turb4 = zeros(jmax,imax);
 
dudy  = zeros(jmax,imax);
muT   = zeros(jmax,imax);
muLam = zeros(jmax,imax);

% Reynolds stresses
uu = zeros(jmax,imax);
vv = zeros(jmax,imax);
ww = zeros(jmax,imax);
uv = zeros(jmax,imax);

% Temporary holder
temp = zeros(imax,size(Data,2));

% Sort into structured (i,j) format
for j = 1:jmax
    temp(1:imax,:) = Data(((j-1)*imax+1):(j*imax),:);
    temp = sortrows(temp,[2]);
    
    for i =1:imax
        XY(j,i,1) = temp(i,1);
        XY(j,i,2) = temp(i,2);
        
        scalars(j,i,1) = temp(i,3);
        scalars(j,i,2) = temp(i,4)./scalars(j,i,1);
        scalars(j,i,3) = temp(i,5)./scalars(j,i,1);
        scalars(j,i,4) = temp(i,6);
        
        if (strcmp(model,'sa'))
            
            turb1(j,i,1) = temp(i,7);
            muLam(j,i,1) = temp(i,8);
            dudy(j,i,1)  = temp(i,9);
            muT(j,i,1)   = temp(i,10);            
            uu(j,i)      = temp(i,11);
            vv(j,i)      = temp(i,12);
            ww(j,i)      = temp(i,13);
            uv(j,i)      = temp(i,14);  
            
        elseif (strcmp(model,'sst') || strcmp(model,'keps'))
            
            turb1(j,i,1) = temp(i,7);
            turb2(j,i,1) = temp(i,8);
            muLam(j,i,1) = temp(i,9);
            dudy(j,i,1)  = temp(i,10);
            muT(j,i,1)   = temp(i,11);
            uu(j,i)      = temp(i,12);
            vv(j,i)      = temp(i,13);
            ww(j,i)      = temp(i,14);
            uv(j,i)      = temp(i,15);
            
        elseif strcmp(model, 'v2f')
            
            turb1(j,i,1) = temp(i,7);
            turb2(j,i,1) = temp(i,8);
            turb3(j,i,1) = temp(i,9);
            turb4(j,i,1) = temp(i,10);
            muLam(j,i,1) = temp(i,11);
            dudy(j,i,1)  = temp(i,12);
            muT(j,i,1)   = temp(i,13);
            uu(j,i)      = temp(i,14);
            vv(j,i)      = temp(i,15);
            ww(j,i)      = temp(i,16);
            uv(j,i)      = temp(i,17);
            
        end
    end
end

clear Data file
end

% Function to decompose matrix of scalars into individual values:
% rho, u, v, P
function [rho,vel,press] = scalar_decomp(scalars)
rho = scalars(:,:,1);
vel(:,:,1) = scalars(:,:,2);
vel(:,:,2) = scalars(:,:,3);
press = scalars(:,:,4);
end

% Function that determines BL thickness
function [delta delta_idx] = blthickness(XY,vel)
delta    = zeros(size(XY,1),1);
delta_idx = zeros(size(XY,1),1);

for i=1:size(XY,1)
    u0 = vel(i,end,1);
    for j=1:size(XY,2)
        if vel(i,j,1) >= .99*u0
            delta(i) = XY(i,j,2);
            delta_idx(i) = j;
            break
        end
    end
    clear u0
end

end

% Function to calculate ReTheta
function [dispth,theta,Re_disp,Re_theta] = calc_Retheta(XY,vel,rho,muLam,freestream,delta_idx)
% For each line of i generate ReTheta based on compressible formulation:
%     Theta = integral( rho(y)/rho0 * u(y)/u0 *(1 - u(y)/u0) dy)

dispth  = zeros(size(XY,1),1);
theta   = zeros(size(XY,1),1);
Re_disp  = zeros(size(XY,1),1);
Re_theta = zeros(size(XY,1),1);

% Calculate for each line i
for i=1:size(XY,1)
    
    y    = XY(i,:,2);    % y values
    u0   = freestream(1); %vel(i,end,1); % free stream velocity
    rho0 = freestream(3); %rho(i,end);   % free stream density
    
    u_u0 = vel(i,:,1) / u0;
    rho_rho0 = rho(i,:) / rho0;
    
    mass = 1-rho_rho0.*u_u0;
    mom = rho_rho0 .* u_u0 .* (1 - u_u0);

    if delta_idx(i)<=1
        dispth(i)  = 0;
        theta(i)   = 0;
        Re_theta(i)= 0;
    else
        dispth(i)  = trapz(y(1:delta_idx(i)),mass(1:delta_idx(i)));
        theta(i)   = trapz(y(1:delta_idx(i)),mom(1:delta_idx(i)));
        Re_disp(i)  = dispth(i,1)*freestream(1)*freestream(3)/freestream(4);
        Re_theta(i) = theta(i,1)*freestream(1)*freestream(3)/freestream(4);
    end
end

end

% Function to extract profile based on some value
function [imatch fak] = extract_point(value)

value_min = min(value);
value_max = max(value);
fprintf('The minimum value is %f\n',value_min);
fprintf('The maximum value is %f\n',value_max);

value_desired = input('For what value would you like to extract?\n','s');
while size(str2num(value_desired)) == 0
    value_desired = input('Input not a real number, try again.\nFor what value value would you like to extract?\n','s');
end
while (str2num(value_desired) < value_min) || (str2num(value_desired) > value_max)
    s = sprintf('Input not in range of data (%f - %f).\nYou may need to re-run the flat-plate simulation',value_min,value_max);
    disp(s);
    value_desired = input('For what value would you like to extract?\n','s');
    while size(str2num(value_desired)) == 0
        value_desired = input('Input not a real number, try again.\nFor what value would you like to extract?\n','s');
    end
end
value_desired = str2num(value_desired);

for i=(length(value)):-1:2
    if value_desired <= value(i) && value_desired > value(i-1)
        imatch = i;
        fak = (value_desired - value(i-1))/(value(i) - value(i-1));
        break;
    end
end
end

% Function to generate the turbulent profile, and output to a file
function gen_profile(XY,vel,rho,press,turb1,turb2,turb3,turb4,muLam,muT,...
    uu,vv,ww,uv,imatch,fak,model)

y(:,1) = XY(imatch-1,:,2) + fak*(XY(imatch,:,2) - XY(imatch-1,:,2));

d(:,1) = rho(imatch-1,:,1) + fak*(rho(imatch,:,1) - rho(imatch-1,:,1));
u(:,1) = vel(imatch-1,:,1) + fak*(vel(imatch,:,1) - vel(imatch-1,:,1));
v(:,1) = vel(imatch-1,:,2) + fak*(vel(imatch,:,2) - vel(imatch-1,:,2));
P(:,1) = press(imatch-1,:) + fak*(press(imatch,:) - press(imatch-1,:));

if strcmp(model,'sa')
    nuhat(:,1) = turb1(imatch-1,:) + fak*(turb1(imatch,:) - turb1(imatch-1,:));
elseif strcmp(model,'sst')
    k(:,1)  = turb1(imatch-1,:) + fak*(turb1(imatch,:) - turb1(imatch-1,:));
    om(:,1) = turb2(imatch-1,:) + fak*(turb2(imatch,:) - turb2(imatch-1,:));
elseif strcmp(model,'v2f')
    k(:,1)   = turb1(imatch-1,:) + fak*(turb1(imatch,:) - turb1(imatch-1,:));
    eps(:,1) = turb2(imatch-1,:) + fak*(turb2(imatch,:) - turb2(imatch-1,:));
    v2(:,1)  = turb3(imatch-1,:) + fak*(turb3(imatch,:) - turb3(imatch-1,:));
    f(:,1)   = turb4(imatch-1,:) + fak*(turb4(imatch,:) - turb4(imatch-1,:));
elseif strcmp(model,'keps')
    k(:,1)   = turb1(imatch-1,:) + fak*(turb1(imatch,:) - turb1(imatch-1,:));
    eps(:,1) = turb2(imatch-1,:) + fak*(turb2(imatch,:) - turb2(imatch-1,:));
end

r11(:,1) = uu(imatch-1,:) + fak*(uu(imatch,:) - uu(imatch-1,:));
r22(:,1) = vv(imatch-1,:) + fak*(vv(imatch,:) - vv(imatch-1,:));
r33(:,1) = ww(imatch-1,:) + fak*(ww(imatch,:) - ww(imatch-1,:));
r12(:,1) = uv(imatch-1,:) + fak*(uv(imatch,:) - uv(imatch-1,:));

eddy_visc(:,1) = muT(imatch-1,:) + fak*(muT(imatch,:) - muT(imatch-1,:));

if (~strcmp(model,'sa'))
    [c1c, c2c, c3c, xbary, ybary] = barymap(r11,r22,r33,r12,d,k);
end

% Variables in wall units
rho_wall   = d(1);
muLam_wall = muLam(imatch-1,1) + fak*(muLam(imatch,1) - muLam(imatch-1,1));
tau_wall   = muLam_wall*(u(1))/(y(1));
u_tau      = sqrt(abs(tau_wall)/rho_wall);
u_plus     = u/u_tau;
y_plus     = (rho_wall*u_tau/muLam_wall)*y;

if (~strcmp(model,'sa'))
    k_plus      = k/(u_tau*u_tau);
end
if strcmp(model,'v2f') || strcmp(model,'keps')
    eps_plus = eps*muLam_wall/(rho_wall*u_tau^4);
end

uprime_plus = zeros(length(y),1);
vprime_plus = zeros(length(y),1);
wprime_plus = zeros(length(y),1);
r12_plus    = zeros(length(y),1);
for i = 1:length(y)
    uprime_plus(i) = sqrt(-1/(d(i)*u_tau*u_tau)*r11(i));
    vprime_plus(i) = sqrt(-1/(d(i)*u_tau*u_tau)*r22(i));
    wprime_plus(i) = sqrt(-1/(d(i)*u_tau*u_tau)*r33(i));
    r12_plus(i)    = -1/(d(i)*u_tau*u_tau)*r12(i);
end
% for i = 1:length(y)
%     uprime_plus(i) = sqrt(r11(i));
%     vprime_plus(i) = sqrt(r22(i));
%     wprime_plus(i) = sqrt(r33(i));
%     r12_plus(i)    = r12(i);
% end

% Output to file:
if strcmp(model,'sa')

    prof = fopen('profilesTec.dat','w');
    
    fprintf(prof,'VARIABLES = "y/h", "u", "v", "P", "nuhat", "y+",');
    fprintf(prof,' "U+", "uprime+", "vprime+", "wprime+", "r12+",');
    fprintf(prof,' "muT"\n');
    
    fprintf(prof,'ZONE t="turbulent profile"\n');
    for jj = 1:length(y)
        j= length(y)-jj+1;
        fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',...
            y(j),u(j),v(j),P(j),nuhat(j),y_plus(j),u_plus(j),...
            uprime_plus(j),vprime_plus(j),wprime_plus(j),r12_plus(j),...
            eddy_visc(i));
    end   
    fclose(prof);
    
    prof = fopen('profilesNew.dat','w');
    fprintf(prof, 'n=%d\td=%d\n', length(y), 6);
    for j = 1:length(y)
        fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',...
            y(j),d(j),u(j),v(j),P(j),nuhat(j));
    end
    fclose(prof);
    
elseif strcmp(model,'sst') 
    
    prof = fopen('profilesTec.dat','w');
    
    fprintf(prof,'VARIABLES = "y/h", "u", "v", "P", "k+", "om", "y+",');
    fprintf(prof,' "U+", "uprime+", "vprime+", "wprime+", "r12+",');
    fprintf(prof,' "muT", "xbary", "ybary"\n');
    
    fprintf(prof,'ZONE t="turbulent profile"\n');
    for jj = 1:length(y)
        j= length(y)-jj+1;
        fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',...
            y(j),u(j),v(j),P(j),k_plus(j),om(j),y_plus(j),u_plus(j),...
            uprime_plus(j),vprime_plus(j),wprime_plus(j),r12_plus(j),...
            eddy_visc(i),xbary(j),ybary(j));
    end   
    fclose(prof);
    
    prof = fopen('profilesNew.dat','w');
    fprintf(prof, 'n=%d\td=%d\n', length(y), 7);
    for j = 1:length(y)
        fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',...
            y(j),d(j),u(j),v(j),P(j),k(j),om(j));
    end
    fclose(prof);
    
elseif strcmp(model,'v2f')
    
    prof = fopen('profilesTec.dat','w');
    
    fprintf(prof,'VARIABLES = "y/h", "u", "v", "P", "k+", "eps+", "v2", "f",');
    fprintf(prof,' "y+", "U+", "uprime+", "vprime+", "wprime+",');
    fprintf(prof,' "r12+", "muT", "xbary", "ybary"\n');
    
    fprintf(prof,'ZONE t="turbulent profile"\n');
    for jj = 1:length(y)
        j= length(y)-jj+1;
        fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',...
            y(j),u(j),v(j),P(j),k_plus(j),eps_plus(j),v2(j),f(j),y_plus(j),u_plus(j),...
            uprime_plus(j),vprime_plus(j),wprime_plus(j),r12_plus(j),...
            eddy_visc(i),xbary(j),ybary(j));
    end  
    fclose(prof);
    
    prof = fopen('profilesNew.dat','w');
    fprintf(prof, 'n=%d\td=%d\n', length(y), 9);
    for j = 1:length(y)
        fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',...
            y(j),d(j),u(j),v(j),P(j),k(j),eps(j),v2(j),f(j));
    end
    fclose(prof);

elseif strcmp(model,'keps')
    
    prof = fopen('profilesTec.dat','w');
    
    fprintf(prof,'VARIABLES = "y/h", "u", "v", "P", "k+", "eps+", "y+",');
    fprintf(prof,' "U+", "uprime+", "vprime+", "wprime+", "r12+",');
    fprintf(prof,' "muT", "xbary", "ybary"\n');
    
    fprintf(prof,'ZONE t="turbulent profile"\n');
    for jj = 1:length(y)
        j= length(y)-jj+1;
        fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',...
            y(j),u(j),v(j),P(j),k_plus(j),eps_plus(j),y_plus(j),u_plus(j),...
            uprime_plus(j),vprime_plus(j),wprime_plus(j),r12_plus(j),...
            eddy_visc(i),xbary(j),ybary(j));
    end  
    fclose(prof);
    
    prof = fopen('profilesNew.dat','w');
    fprintf(prof, 'n=%d\td=%d\n', length(y), 7);
    for j = 1:length(y)
        fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',...
            y(j),d(j),u(j),v(j),P(j),k(j),eps(j));
    end
    fclose(prof);    
end

end

% Function to generate wall values
function [Re_tau] = gen_wall_values(wallfile,freestream,Re_theta,delta)
    
% The following comes from wall_values.m    
file = wallfile;
raw = importdata(file);
wall_val = sortrows(raw.data);

Uinf   = freestream(1);
P01    = freestream(2);
rhoinf = freestream(3);

x        = wall_val(:,1);
rho      = wall_val(:,4);
P        = wall_val(:,5);
tau_wall = wall_val(:,7);
yplus    = wall_val(:,9);
mulam    = wall_val(:,10);

utau  = sqrt(abs(tau_wall)./rho);
Cf    = tau_wall/(0.5*rhoinf*Uinf^2);
Cp    = (P-P01)/(0.5*rhoinf*Uinf^2);

diff = length(Re_theta)-length(x);
wallRe_theta = Re_theta(diff+1:end);

wallRe_tau = zeros(length(utau),1);
for i = 1:length(utau)
    wallRe_tau(i) = rhoinf*utau(i)*delta(i+diff)/mulam(i);
end
Re_tau = zeros(1,length(Re_theta));
Re_tau(diff+1:end) = wallRe_tau;


% Output to file: x Cf Cp yplus Retau wallRetheta
prof = fopen('wall_values.dat','w');
fprintf(prof,'VARIABLES = "x", "Cf", "Cp", "yplus", "Retau", "wallRetheta"\n');
fprintf(prof,'ZONE t="wall values"\n');
for j = 1:length(x)
    fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',...
        x(j),Cf(j),Cp(j),yplus(j),wallRe_tau(j),wallRe_theta(j));
end
fclose(prof);
end

% Function to generate barycentric map
function [c1c, c2c, c3c, xbary, ybary] = barymap(R11,R22,R33,R12,d,k)
a = zeros(3,3);
a11(:,1) = -R11./(2*d.*k) - 1.0/3.0;
a22(:,1) = -R22./(2*d.*k) - 1.0/3.0;
a33(:,1) = -R33./(2*d.*k) - 1.0/3.0;
a12(:,1) = -R12./(2*d.*k);

c1c = zeros(length(d),1);
c2c = zeros(length(d),1);
c3c = zeros(length(d),1);

xbary = zeros(length(d),1);
ybary = zeros(length(d),1);

for i = 1:length(d)
    a(1,1) = a11(i);
    a(1,2) = a12(i);
    a(1,3) = 0.0;
    
    a(2,1) = a(1,2);
    a(2,2) = a22(i);
    a(2,3) = 0.0;
    
    a(3,1) = a(1,3);
    a(3,2) = a(2,3);
    a(3,3) = a33(i);
    
    lambda = sort(eig(a),'descend');
    
    c1c(i) = lambda(1) - lambda(2);
    c2c(i) = 2.0*(lambda(2) - lambda(3));
    c3c(i) = 3.0*lambda(3) + 1.0;
    
    xbary(i) = c1c(i)*1.0 + c2c(i)*0.0 + c3c(i)*0.5;
    ybary(i) = c1c(i)*0.0 + c2c(i)*0.0 + c3c(i)*0.866025;
end
end


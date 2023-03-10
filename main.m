clear
clc
close all

%% Global parameters
mu = 39.4769; %Sun's gravitational parameter, (au^3/year^2)

%unit conversions
meter = 1/149597870700; %1m to au relation
second = 1/(86400*365.25); %1 second to year relation

%initial coordinates of spaceship (Sun-centered intertial frame)
r0 = [-1.05;0;0]; %au
v0 = [0;-6.1316;0]; %au/year ^j

%acceleration from propulsion system
aT0 = (1/3) * 1e-4; % m*s^-2
aT0 = aT0 * meter/(second^2);%au/year

ad_vect =  @(r_mag,v_unit) aT0 .* (1./r_mag).^2 .* (v_unit);


%% Part 1 - Cowell's Method
%%%
% Directly solving perturbed two-body problem with ode45
%%%

% r_dd + mu*r_vect/(r_mag.^3)




tspan = [0,20];
y0 = [r0;v0];
[t,r] = ode45(@(t,r) cowell(r,mu,ad_vect),tspan,y0);

figure()
plot(r(:,1),r(:,2),'-x')
hold on
plot(0,0,'ro') %the sun



%% Part 2 - Encke's Method


%% Part 3 - Osculating Elements


%% Function definitions
% Part 1
function stateSpaceRepCowell = cowell(y,mu,ad_vect_fun)
    r = y(1:3);
    v = y(4:6);
    
    v_mag = sum(v.^2);
    v_unit = v./v_mag;
%     disp("Unit v")
%     disp(v_unit)
    r_mag = sum(r.^2);
%     r_unit = r(1)/r_mag;
    ad_vect = ad_vect_fun(r_mag,v_unit);
%     disp("Perturbation vector");
%     disp(ad_vect)
    stateSpaceRepCowell = [v;...
                           -(mu./(r_mag.^3)).*r + ad_vect];
end
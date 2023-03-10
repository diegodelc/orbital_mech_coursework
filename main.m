clear
clc
close all

%% Global parameters
mu = 39.4769; %Sun's gravitational parameter, (au^3/year^2)

%unit conversions
one_meter = 1/149597870700; %1m to au relation
one_second = 1/(86400*365.25); %1 second in years
acc_to_au_and_years = 6656.77641; %1m/s^2 to au/year^2


%initial coordinates of spaceship (Sun-centered intertial frame)
r0 = [-1.05;0;0]; %au
v0 = [0;-6.1316;0]; %au/year ^j

%acceleration from propulsion system
aT0 = (1/3) * 10^-4; %m*s^-2
aT0 = aT0 * one_meter/(one_second^2); %au/year^2
% aT0 = aT0 * acc_to_au_and_years;


ad_vect =  @(r_mag,v_unit) aT0 * ((1./r_mag).^2 ).* (v_unit);


%% Part 1 - Cowell's Method
%%%
% Directly solving perturbed two-body problem with ode45
%%%



tspan = [0,20];
y0 = [r0;v0];
[t,y] = ode45(@(t,y) cowell(y,mu,ad_vect),tspan,y0);



val = 0;
figure()
plot(y(:,1+val),y(:,2+val),'-.')
hold on
plot(0,0,'ro') %the sun
axis equal
xlabel("i")
ylabel("j")


%% Part 2 - Encke's Method


%% Part 3 - Osculating Elements


%% Function definitions
% Part 1
function stateSpaceRepCowell = cowell(y,mu,ad_fun)
    r = y(1:3);
    v = y(4:6);
    
    
    v_mag = sum(v.^2);
    v_unit = v./v_mag;

    r_mag = sum(r.^2);
%     r_unit = r./r_mag;
    
    ad = ad_fun(r_mag,v_unit);

    stateSpaceRepCowell = [v;ad-(mu.*r)./(r_mag.^3)];
end
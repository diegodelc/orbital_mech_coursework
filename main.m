clear
clc
close all

%% Global parameters
mhu = 39.4769; %Sun's gravitational parameter, (au^3/year^2)



%initial coordinates of spaceship (Sun-centered intertial frame)
r0 = -1.05; %au î (i unit vector)
v0 = -6.1316; %au/year ^j (j unit vector)

%acceleration from propulsion system
aT0 = (1/3) * 1e-4;
aT_vect =  @(r,v_unit) aT0 .* (1/r).^2 .* (v_unit);


%% Part 1 - Cowell's Method


%% Part 2 - Encke's Method


%% Part 3 - Osculating Elements
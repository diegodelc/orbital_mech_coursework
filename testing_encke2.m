clear
clc
close all

%% Global parameters
global mu
mu = 39.4769; %Sun's gravitational parameter, (au^3/year^2)

%unit conversions
one_meter = 1/149597870700; %1m to au relation
one_second = 1/(86400*365.25); %1 second in years
acc_to_au_and_years = 6656.77641; %1m/s^2 to au/year^2
hours = 3600; %Hours to seconds
days = 24*hours; %Days to seconds
deg = pi/180;


%initial coordinates of spaceship (Sun-centered inertial frame)
r0 = [-1.05;0;0]; %au
v0 = [0;-6.1316;0]; %au/year ^j

%acceleration from propulsion system
aT0 = (1/3) * 10^-4; %m*s^-2
aT0 = aT0 * one_meter/(one_second^2); %au/year^2
% aT0 = aT0 * acc_to_au_and_years;

%function handle for calculating acceleration
ad_vect =  @(r_mag,v_unit) aT0 * ((1./r_mag).^2 ).* (v_unit);





%% Part 1 - Cowell's Method
%%%
% Directly solving perturbed two-body problem with ode45
%%%
del_t = 0.01;
options = odeset('maxstep', del_t);
tspan = [0,20];
y0 = [r0;v0];
[t_cowell,y_cowell] = ode45(@(t,y) cowell(y,ad_vect),tspan,y0,options);

%% Part 2 - Encke's Method
%%%
% Solving perturbed two-body problem using Encje's method
%   Code implementation is based on framework from [1]

% [1] Orbital Mechanics for Engineering Students, Fourth Edition, Howard D.Curtis, 
%%%


t = [0,20]; %in years
t0 = t(1);
tf = t(2);

R0 = r0';

r0 = norm(R0);

V0 = v0';
v0 = norm(V0);

%Time step for Encke procedure
del_t = 5/365; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = odeset('maxstep', del_t);

%Begin Encke integration
t = t0;
tsave = t0;
y = [R0 V0];

del_y0 = zeros(1,6);

t = t + del_t;

while t <= tf + del_t/2

    [time,z] = ode45(@(t,f) rates_diego(t,f,R0,V0,t0,ad_vect,del_t), [t0 t], del_y0, options);


    [Rosc,Vosc] = propDt(R0, V0, t,t0,del_t);

    
    R0 = Rosc + z(end,1:3);
    
    V0 = Vosc + z(end,4:6);
    


    t0 = t;
    tsave = [tsave;t];
    y = [y; [R0 V0]];
    t = t + del_t;
    del_y0 = zeros(1,6);

end
y_encke_diego = y;


%% Plotting
% % Plotting Cowell's method only
% val = 0;
% figure()
% plot(0,0,'ro','DisplayName', 'SUN') %the sun
% hold on
% plot(y_cowell(:,1+val),y_cowell(:,2+val),'-.b','DisplayName',"Cowell's Method")
% legend
% axis equal
% title("Cowell's method")
% 
% %Plot Encke's method only
% figure()
% plot(0,0,'ro','DisplayName', 'SUN')
% hold on
% plot(y_encke(:,1),y_encke(:,2),'-.b','DisplayName', "Encke's Method")
% legend
% axis equal
% title("Encke's Methd")

%Plotting both methods together
figure()
plot(0,0,'ro','DisplayName', 'SUN')
hold on
plot(y_cowell(:,1),y_cowell(:,2),'-.b','DisplayName',"Cowell's Method")
hold on
plot(y_encke_diego(:,1),y_encke_diego(:,2),'-kx','DisplayName', "Encke's Method")
legend
axis equal
title("Cowell's and Encke's Methds")

%% Part 3 - Osculating Elements
t = tsave;
n_times = length(t);

for j = 1:n_times
    R = [y_encke_diego(j,1:3)];
    V = [y_encke_diego(j,4:6)];
    r(j) = norm(R);
    v(j) = norm(V);
    coe = coe_from_sv(R,V, mu);
    h(j) = coe(1);
    e(j) = coe(2);
    RA(j) = coe(3);
    i(j) = coe(4);
    w(j) = coe(5);
    TA(j) = coe(6);
end

figure()

subplot(2,1,1)
plot(t/3600,(RA)/deg)
title('Variation of Right Ascension')
xlabel('hours')
ylabel('{\it\Delta\Omega} (deg)')
grid on
grid minor
axis tight

subplot(2,1,2)
plot(t/3600,(w)/deg)
title('Variation of Argument of Perigee')
xlabel('hours')
ylabel('{\it\Delta\omega} (deg)')
grid on
grid minor
axis tight

figure()
subplot(3,1,1)
plot(t/3600,h)
title('Variation of Angular Momentum')
xlabel('hours')
ylabel('{\it\Deltah} (km^2/s)')
grid on
grid minor
axis tight

subplot(3,1,2)
plot(t/3600,e)
title('Variation of Eccentricity')
xlabel('hours')
ylabel('\it\Deltae')
grid on
grid minor
axis tight

subplot(3,1,3)
plot(t/3600,(i)/deg)
title('Variation of Inclination')
xlabel('hours')
ylabel('{\it\Deltai} (deg)')
grid on
grid minor
axis tight

%% Function definitions
% Part 1: Cowell's Method (function(s))
function stateSpaceRepCowell = cowell(y,ad_fun)
    r = y(1:3);
    v = y(4:6);
    global mu
    
    v_mag = norm(v);
    v_unit = v/v_mag;

    r_mag = norm(r);
%     r_unit = r./r_mag;
    
    ad = ad_fun(r_mag,v_unit);

    stateSpaceRepCowell = [v;ad-(mu.*r)/(r_mag^3)];
end

%Part 2: Encke's method (functions)
function coe = coe_from_sv(R,V,mu)
    eps = 1.e-10;
    
    r = norm(R);
    v = norm(V);

    vr = dot(R,V)/r;
    H = cross(R,V);
    h = norm(H);
    
    incl = acos(H(3)/h);
    
    N = cross([0 0 1],H);
    n = norm(N);
    
    if n ~= 0
        RA = acos(N(1)/n);
        if N(2) < 0
            RA = 2*pi - RA;
        end
    else
    RA = 0;
    end
    
    E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
    e = norm(E);
    
    if n ~= 0
        if e > eps
            w = acos(dot(N,E)/n/e);
            if E(3) < 0
                w = 2*pi - w;
            end
        else
            w = 0;
        end
    else
        w = 0;
    end
    
    if e > eps
        TA = acos(dot(E,R)/e/r);
        if vr < 0
            TA = 2*pi - TA;
        end
    else
        cp = cross(N,R);
        if cp(3) >= 0
            TA = acos(dot(N,R)/n/r);
        else
            TA = 2*pi - acos(dot(N,R)/n/r);
        end
    end
    
    a = h^2/mu/(1 - e^2);
    coe = [h e RA incl w TA a];
end

function dfdt = rates_diego(t,f,R0,V0,t0,ad_vect,del_t)
    del_r = f(1:3)'; %Position deviation
    del_v = f(4:6)'; %Velocity deviation
    
%     disp("from rates: ")
    [Rosc,Vosc] = propDt(R0, V0, t,t0,del_t);
    
   
        
    Rpp = Rosc + del_r;
    Vpp = Vosc + del_v;
%     Rpp = Rosc;
%     Vpp = Vosc;
%     rosc = norm(Rosc);
%     rpp = norm(Rpp);
    
    %compute perturbing acceleration (here for J2, need to change to ours)
%     xx = Rpp(1); yy = Rpp(2); zz = Rpp(3);
%     fac = 3/2*J2*(mu/rpp^2)*(RE/rpp)^2;
%     ap = -fac*[(1 - 5*(zz/rpp)^2)*(xx/rpp) ...
%         (1 - 5*(zz/rpp)^2)*(yy/rpp) ...
%         (3 - 5*(zz/rpp)^2)*(zz/rpp)];
% 
%     F = 1 - (rosc/rpp)^3;
%     del_a = -mu/rosc^3*(del_r - F*Rpp) + ap;
    v_in = Vpp;
    v_unit = v_in/norm(v_in);
    
    r_in = Rpp;
    r_mag = norm(r_in);
    del_a = ad_vect(r_mag,v_unit);
    %are these two the same?
%     dfdt = [del_v(1) del_v(2) del_v(3) del_a(1) del_a(2) del_a(3)]';
    dfdt = [del_v del_a]';

end 


function [R,V] = propDt(R0,V0,t,t0,del_t)
    global mu

    r0 = norm(R0);
    v0 = norm(V0);
    
    %semi-major axis (a)
    a = -mu/(v0^2-2*mu/r0);

    %find angular momentum (h)
    H = cross(R0,V0);
    h = norm(H);

    %find eccentricity (e)
    E = (1/mu)*cross(V0,H) - R0/r0;
    e = norm(E);

    %find theta0
    theta0 = acos(dot(R0,E)/(r0*e));
    
    
    %find E0
    E0 = 2*atan(sqrt((1-e)/(1+e))*tan(theta0/2));


    %find M0
    M0 = E0 - e*sin(E0);

    %find n
    n = sqrt(mu/a^3);

    %find M1
    M1 = M0 + n*del_t;

    %iterate to find E
    E1 = M1;
    mydiff = 1; %dummy value to get into while loop
    counter = 0;
    while abs(mydiff) > 1e-8
%     for i = 1:5
        mydiff = E1;
        E1 = M1 + e*sin(E1);
        mydiff = mydiff-E1;
        counter = counter +1;
    end
%     counter
    
    %find theta1
    theta1 = 2*atan(sqrt((1+e)/(1-e))*tan(E1/2)); %if e < 1 -> imaginary number!!
    del_theta = theta1-theta0;

    %find r1
    P = h^2/mu;
    r = P/(1+e*cos(theta1));

    %find f and g
    f = 1 - (mu*r/h^2)*(1-cos(del_theta));
    g = (r*r0/h)*sin(del_theta);

    %find fdot and gdot
    gdot = 1 - (mu*r0/h^2)*(1-cos(del_theta));
%     fdot = (mu/h) * ...
%         ((1-cos(del_theta))/sin(del_theta)) * ...
%         ((mu/h^2)*(1-cos(del_theta)) - (1/r0) - (1/r1));
    fdot = (f*gdot-1)/g;

    %find R and V
    R = f*R0 + g*V0;
    V = fdot*R0 + gdot*V0;

    
    
%     R = R';
%     V = V';

end


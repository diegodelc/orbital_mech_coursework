clear
close all


%% Global
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


ad_vect =  @(r_mag,v_unit) aT0 * ((1./r_mag).^2 ).* (v_unit);


%% Encke's Method



t = [0,20]; %in years
t0 = t(1);
tf = t(2);

R0 = r0';

r0 = norm(R0);

V0 = v0';
v0 = norm(V0);

%Time step for Encke procedure
del_t = 0.1;
options = odeset('maxstep', del_t);

%Begin Encke integration
t = t0;
tsave = t0;
y = [R0 V0];

del_y0 = zeros(1,6);

t = t + del_t;

while t <= tf + del_t/2
    [dum,z] = ode45(@(t,f) rates(t,f,R0,V0,t0,ad_vect), [t0 t], del_y0, options);

    [Rosc,Vosc] = rv_from_r0v0(R0', V0', t-t0);

    R0 = Rosc' + z(end,1:3);
    
    V0 = Vosc' + z(end,4:6);

    t0 = t;
    tsave = [tsave;t];
    y = [y; [R0 V0]];
    t = t + del_t;
    del_y0 = zeros(1,6);
end
%% Plotting osculating elements
t = tsave;
n_times = length(t);

for j = 1:n_times
    R = [y(j,1:3)];
    V = [y(j,4:6)];
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

%% Plotting orbit
figure(1)
plot(0,0,'ro','DisplayName', 'SUN')
hold on
plot(y(:,1),y(:,2),'--bx','DisplayName', "Encke's Method")
legend
axis equal

%%
figure(2)

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

figure(3)
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
%%


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
    %...Equation 4.13a (incorporating the case e = 0):
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

function dfdt = rates(t,f,R0,V0,t0,ad_vect)
    del_r = f(1:3)'; %Position deviation
    del_v = f(4:6)'; %Velocity deviation
    
    
    [Rosc,Vosc] = rv_from_r0v0(R0, V0, t-t0);
    
    Rpp = Rosc + del_r;
    Vpp = Vosc + del_v;
%     rosc = norm(Rosc);
%     rpp = norm(Rpp);
    
    %compute perturbing acceleration (here for J2, need to change to ours)

%     RE = earthRadius;
%     xx = Rpp(1); yy = Rpp(2); zz = Rpp(3);
%     fac = 3/2*J2*(mu/rpp^2)*(RE/rpp)^2;
%     ap = -fac*[(1 - 5*(zz/rpp)^2)*(xx/rpp) ...
%         (1 - 5*(zz/rpp)^2)*(yy/rpp) ...
%         (3 - 5*(zz/rpp)^2)*(zz/rpp)];
% 
%     F = 1 - (rosc/rpp)^3;
%     del_a = -mu/rosc^3*(del_r - F*Rpp) + ap;
    
    del_a = ad_vect(norm(Rpp),Vpp/norm(Vpp));
    %are these two the same?
    dfdt = [del_v(1) del_v(2) del_v(3) del_a(1) del_a(2) del_a(3)]';
%     dfdt = [del_v del_a]';

end 

function s = stumpS(z)
    if z>0
        s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif z<0
        s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    else
        s = 1/6;
    end
end

function c = stumpC(z)
    if z > 0
        c = (1 - cos(sqrt(z)))/z;
    elseif z < 0
        c = (cosh(sqrt(-z)) - 1)/(-z);
        %%e44 MATLAB scripts
    else
        c = 1/2;
    end
end

function x = kepler_U(dt, ro, vro, a)
    global mu
    
    error = 1e-8;
    nMax = 1000;

    
    x = sqrt(mu)*abs(a)*dt;
    
    n=0;
    ratio = 1;
    while abs(ratio) > error && n <= nMax
        
        n= n+1;
        C = stumpC(a*x^2);
        S = stumpS(a*x^2);

        
        F = ro*vro/sqrt(mu)*x^2*C + (1 - a*ro)*x^3*S + ro*x - sqrt(mu)*dt;
        dFdx = ro*vro/sqrt(mu)*x*(1 - a*x^2*S) + (1 - a*ro)*x^2*C + ro;
        ratio = F/dFdx;
        x = x - ratio;
    end
    
    if n > nMax
        fprintf("\n **No. iterations of Kepler's equation = %g", n)
        fprintf('\n F/dFdx = %g\n', F/dFdx)
    end
    
end

function [f,g] = f_and_g(x,t,ro,a)
    global mu
    
    z = a*x^2;

    f = 1 - x^2/ro*stumpC(z);
    
    g = t - 1/sqrt(mu)*x^3*stumpS(z);
    
    
end

function [fdot, gdot] = fDot_and_gDot(x, r, ro, a)
    global mu
    z = a*x^2;
    
    fdot = sqrt(mu)/r/ro*(z*stumpS(z) - 1)*x;
    
    gdot = 1 - x^2/r*stumpC(z);
end

function [R,V] = rv_from_r0v0(R0,V0,t)
    global mu
    
    r0 = norm(R0);
    v0 = norm(V0);
    
    %initial radial velocity
    vr0 = dot(R0,V0)/r0;
    
    %Reciprocal of the semimajor axis
    alpha = 2/r0 - v0^2/mu;
    
    %Universal anomaly
    
    x = kepler_U(t,r0,vr0,alpha);
    
    %f and g functions:
    [f, g] = f_and_g(x, t, r0, alpha);
    
    %Final position vector
    R = f*R0 + g*V0;
    
    %Magnitude of R
    r = norm(R);
    
    %Derivatives of f and g
    [fdot, gdot] = fDot_and_gDot(x, r, r0, alpha);
    
    %...Compute the final velocity:
    V = fdot*R0 + gdot*V0;
end


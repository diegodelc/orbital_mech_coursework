clear
clc
close all

%% Global parameters
global mu
mu = 39.4769; %Sun's gravitational parameter, (au^3/year^2)

%unit conversions
one_meter = 1/149597870700; %1m to au
one_second = 1/(86400*365.25); %1 second to years
hours = 3600; %Hours to seconds
days = 24*hours; %Days to seconds
deg = pi/180;


%initial coordinates of spaceship (Sun-centered inertial frame)
R0_og = [-1.05;0;0]; %au
V0_og = [0;-6.1316;0]; %au/year ^j

%acceleration from propulsion system
aT0 = (1/3) * 10^-4; %m*s^-2
aT0 = aT0 * one_meter/(one_second^2); %au/year^2


%function handle for calculating acceleration
ad_vect =  @(r_mag,v_unit) aT0 * ((1./r_mag).^2 ).* (v_unit);

tspan = [0,20]; %years

%% Part 1 - Cowell's Method
%%%
% Directly solving perturbed two-body problem with ode45
%%%

del_t = 1/365.25; %years
options = odeset('maxstep', del_t);
% tspan = [0,20];
y0 = [R0_og;V0_og];
[t_cowell,y_cowell] = ode45(@(t,y) cowell(y,ad_vect),tspan,y0,options);

final = y_cowell(end,:);
rfinal = final(1:3)
vfinal = final(4:6)

%% Part 2 - Encke's Method (using eccentric anomaly)
%Time step for both Encke procedures
del_t = 20; %time is in years
options = odeset('maxstep', del_t,'abstol',1e-8,'reltol',1e-8);


%These will be overwritten
t0 = tspan(1);
tf = tspan(2);

R0 = R0_og';
r0 = norm(R0);

V0 = V0_og';
v0 = norm(V0);


%Begin Encke integration
t = t0;
tsave = t0;
y = [R0 V0];

del_y0 = zeros(1,6);

t = t + del_t;
saved_states = zeros(1,6);
while t < tf
    
    [time,z] = ode45(@(t,f) rates_diego(t,f,R0,V0,t0,ad_vect,del_t), [t0 t], del_y0, options);


    [Rosc,Vosc] = propDt(R0, V0,del_t);

    
    R0 = Rosc + z(end,1:3);
    
    V0 = Vosc + z(end,4:6);
    


    t0 = t;
    tsave = [tsave;t];
    y = [y; [R0 V0]];
    
    t = t + del_t;
    del_y0 = zeros(1,6);
    
   

end

if t0 ~= tf
    t = tf;
    [time,z] = ode45(@(t,f) rates_diego(t,f,R0,V0,t0,ad_vect,del_t), [t0 t], del_y0, options);
    
    
    [Rosc,Vosc] = propDt(R0, V0,del_t);
    
    
    R0 = Rosc + z(end,1:3);
    
    V0 = Vosc + z(end,4:6);
    
    
    
%     t0 = t;
    tsave = [tsave;t];
    y = [y; [R0 V0]];
    
%     t = t + del_t;
%     del_y0 = zeros(1,6);
end
tsave_encke_diego = tsave;
y_encke_diego = y;


%% Part 2 - Encke's Method (using unversal anomaly)
%%%
% Solving perturbed two-body problem using Encke's method
%   Code implementation is based on framework from [1]

% [1] Orbital Mechanics for Engineering Students, Fourth Edition, Howard D.Curtis, 
%%%



%These will be overwritten
t0 = tspan(1);
tf = tspan(2);

R0 = R0_og';
r0 = norm(R0);

V0 = V0_og';
v0 = norm(V0);



%Begin Encke integration
t = t0;
tsave = t0;
y = [R0 V0];

del_y0 = zeros(1,6);

t = t + del_t;

while t <= tf + del_t/2
    [~,z] = ode45(@(t,f) rates_book(t,f,R0,V0,t0,ad_vect), [t0 t], del_y0, options);

    [Rosc,Vosc] = rv_from_r0v0(R0, V0, t-t0);

    R0 = Rosc + z(end,1:3);
    
    V0 = Vosc + z(end,4:6);

    t0 = t;
    tsave = [tsave;t];
    y = [y; [R0 V0]];
    t = t + del_t;
    del_y0 = zeros(1,6);
end
tsave_encke = tsave;
y_encke = y;


%% Plotting Orbits

% %Plotting Cowell's method
% figure()
% plot(0,0,'ro','DisplayName', 'SUN')
% hold on
% plot(y_cowell(:,1),y_cowell(:,2),'-.b','DisplayName',"Cowell's Method")
% legend
% xlabel("$\hat{i}$ \textit{(AU)}","interpreter","latex")
% ylabel("$\hat{j}$ \textit{(AU)}","interpreter","latex")
% axis equal
% % title("Cowell's Method")

%Plotting all methods together
figure()
plot(0,0,'go','DisplayName', 'SUN')
hold on
plot(y_cowell(:,1),y_cowell(:,2),'-.b','DisplayName',"Cowell's Method")
hold on
plot(y_encke(:,1),y_encke(:,2),'-xk','DisplayName', "Encke's Method - universal anomaly")
hold on
plot(y_encke_diego(:,1),y_encke_diego(:,2),'--or','DisplayName', "Encke's Method - eccentric anomaly")
xlabel("$\hat{i}$ \textit{(AU)}","interpreter","latex")
ylabel("$\hat{j}$ \textit{(AU)}","interpreter","latex")
legend
axis equal
% title("Cowell's and Encke's Methods")

%% Part 3 - Osculating Elements


times = t_cowell(1:length(t_cowell-1)/4:end)
times = [times;t_cowell(end)]
ys = [R0_og' V0_og'];

for i = 2:length(times)
    ys = [ys; y_cowell(times(i) == t_cowell,:)];
    
end

ys;
[row,col] = size(ys);
figure()
plot(0,0,'go')
hold on
for i = 1:row
    
    R0 = ys(i,1:3);
    V0 = ys(i,4:6);
    sol = [R0 V0];
    for j = 1:5/365:10
        
        [R0,V0] = propDt(R0,V0,5/365); %% propagate each five years
        sol = [sol;[R0 V0]];
    end
    plot(sol(:,1),sol(:,2),'-','DisplayName', sprintf("Initial time: %.1f years",times(i)))
    hold on
end
plot(y_cowell(:,1),y_cowell(:,2),'-c','DisplayName',"Cowell's Method")
hold on

%plot the starting points
for i=1:5
    plot(ys(i,1),ys(i,2),'x','DisplayName',"")
    hold on
end
legend("Sun",...
       sprintf("Initial time: %.1f years",times(1)),...
       sprintf("Initial time: %.1f years",times(2)),...
       sprintf("Initial time: %.1f years",times(3)),...
       sprintf("Initial time: %.1f years",times(4)),...
       sprintf("Initial time: %.1f years",times(5)),...
       "Cowell's Method",...
       'Location', 'best');
axis equal
xlabel("$\hat{i}$ \textit{(AU)}","interpreter","latex")
ylabel("$\hat{j}$ \textit{(AU)}","interpreter","latex")
sol;
%%





which_sim = y_encke_diego;
sim_times = tsave_encke_diego;
% which_sim = y_encke;
% which_sim = y_cowell;

[tsteps,~] = size(which_sim);
a = zeros(tsteps,1);
for j = 1:tsteps
    R = [which_sim(j,1:3)];
    V = [which_sim(j,4:6)];
%     r(j) = norm(R);
%     v(j) = norm(V);
    coe = coe_from_sv(R,V, mu);
%     h(j) = coe(1);
    e(j) = coe(2);
%     RA(j) = coe(3);
%     i(j) = coe(4);
%     w(j) = coe(5);
%     TA(j) = coe(6);
    a(j) = coe(7);
end


%same plots but picking only 5 pointsç
n = 5; %five points
step = floor(tsteps/(n-1));

figure()
% plot(sim_times(1:step:end),a(1:step:end)','--x')
% hold on
% myFit = fit(sim_times(1:step:end),a(1:step:end),'poly1');
% plot(myFit)
% hold on
plot(sim_times,a,'-.k') %all data in cyan (very light)
% title("Variation of semi-major axis")
xlabel('\itTime (years)')
ylabel('\itsemi-major axis (AU)')
% legend("semi-major axis","linear fit","Location","NorthWest")

figure()
% plot(sim_times(1:step:end),e(1:step:end),'--x')
% hold on
% myFit = fit(sim_times(1:step:end),e(1:step:end)','poly1');
% plot(myFit)
% hold on
plot(sim_times,e,'-.k') %all data in cyan (very light)
% title('Variation of Eccentricity')
xlabel('\itTime (years)')
ylabel('\it\Deltae')
% legend("eccentricity","linear fit", "Location", "NorthWest")
% 

%% Function definitions
% Part 1: Cowell's Method (function)
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
function dfdt = rates_book(t,f,R0,V0,t0,ad_vect)
    global mu
    del_r = f(1:3)'; %Position deviation
    del_v = f(4:6)'; %Velocity deviation
    
    
    [Rosc,Vosc] = rv_from_r0v0(R0, V0, t-t0);
    
    Rpp = Rosc + del_r;
    Vpp = Vosc + del_v;

    v_in = Vpp;
    v_unit = v_in/norm(v_in);
    
    r_in = Rpp;
    r_mag = norm(r_in);
    
    
    del_a = (-mu/norm(Rosc)^3)*(del_r- (1-norm(Rosc)^3/norm(Rpp)^3).*Rpp) + ... 
            ad_vect(r_mag,v_unit);

    dfdt = [del_v del_a]';

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
    
    %Compute the final velocity:
    V = fdot*R0 + gdot*V0;
end

% Part 2 - my way
function dfdt = rates_diego(t,f,R0,V0,t0,ad_vect,del_t)
    global mu
    del_r = f(1:3)'; %Position deviation
    del_v = f(4:6)'; %Velocity deviation
    
%     disp("from rates: ")
    [Rosc,Vosc] = propDt(R0, V0,t-t0);
    
   
        
    Rpp = Rosc + del_r;
    Vpp = Vosc + del_v;

%     v_in = Vpp;
%     v_unit = v_in/norm(v_in);
    
%     r_in = Rpp;
%     r_mag = norm(r_in);
    
    ad = ad_vect(norm(Rpp),Vpp/norm(Vpp));
    q = dot(del_r, del_r - 2*Rpp)/(norm(Rpp)^2);
    F = -q*(3+3*q+q^2) / (1+(1+q)^(1.5));
    del_a = -mu*del_r/(norm(Rosc)^3) + mu*Rpp/(norm(Rosc)^3)*F + ad;


    dfdt = [del_v del_a]';

end 


function [R,V] = propDt(R0,V0,del_t)
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
    if abs(del_theta) < 1e-8
       fdot = 0;
    else
        fdot = (mu/h) * ...
            ((1-cos(del_theta))/sin(del_theta)) * ...
            ((mu/h^2)*(1-cos(del_theta)) - (1/r0) - (1/r));
    %     fdot = (f*gdot-1)/g;
    end
    %find R and V
    R = f*R0 + g*V0;
    V = fdot*R0 + gdot*V0;

    
    
%     R = R';
%     V = V';

end

%Part 3 - Osculating elements (function)
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

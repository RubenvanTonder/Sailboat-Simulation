%%% Numerical Parameter Calculation for dragonflite95
% Water Density
rho = 1025;

%% Ship Dimensions
B = 0.125;  % Breadth
T = 0.4;    % Draft
L = 0.95;   % Length
m = 2.7;      % Mass
CB = 0.05;  % Centre of Buoyancy
U = 1.3;    % Hull Speed
u = 1.3;    % Boat Speed

a = L/2;    
b = B/2;    
c = T;  

e = sqrt(1-(b^2)/(a^2));
Ao = (2*(1-e^2))/(e^3) * ((1/2) * log((1+e)/(1-e))-e);
Bo = 1/e^2 - (1-e^2)/(2*e^2)*log((1+e)/(1-e));
Co = Bo;

% Displaced Water Volume
Volume = 4/3 * pi * a * b *c /2;

%% Hydrodynmaic Coefficients
k11 = (Ao)/(2-Ao);
k22 = (Bo)/(2-Bo);
k33 = (Co)/(2-Co);
k44 = 0;
k55 = ((L^2-4*T^2)^2 *(Ao-Co))/(2*(c^4 - a^4) + (Co-Ao)*(4*T^2+L^2)^2);
k66 = ((L^2-4*B^2)^2 *(Bo-Ao))/(2*(L^4-B^4)+(Ao-Bo)*(L^2+B^2)^2);

%% Moment of Inertia of the displaced water
Ixx = 1/120 * pi * rho * L * B * T * (4*T^2 + B^2);
Iyy = 1/120 * pi * rho * L * B * T * (4*T^2 + L^2);
Izz = 1/120 * pi * rho * L * B * T * (B^2 + L^2);

%% Hull Shape Parameters
H = 0.5;
beta = 0.785;

q = (H-1)/(H+1);
p = beta - pi/4;

b = (3/4 * pi *sqrt((pi/4)^2 - pi/2 * p * (1-q^2)))/(pi + p*(1-q^2));
a = (b+1)*q;

% Centre of buoyancy
xb = -0.05;

%% Added Mass Components
m11 = m*k11;
m22 = m*k22;
m44 = k44*Ixx;
m66 = k66*Izz;
m26 = xb*m22;
m62 = m26;
m24 = 0;
m42 = m24;

%% Linear Damping Forces and Moments
normTerm = 1/2 * rho * U^2 * L^2;
Yv = -pi * T^2/L^2 * (1 +0.4 * CB * B/T);
Yr = pi *  T^2/L^2 * (0.5 -2.2 * B/L +0.08* B/T);
Nv =  -pi * T^2/L^2 * (0.5 +2.4 * T/L);
Nr = -pi * T^2/L^2 * (0.25 +0.039 * B/T - 0.56 * B/L);


% Xu damping
Ct = 0.0036 + (0.0152)/(log(L) + 0.6);
Xu = 2*Ct* rho * Volume ^(2/3);


%% Hull Resistance Calculation 
% Reynolds Number
L = 0.95; %m
V = 1.3; %m/sec
v = 3.32926*10^(-6); %m^2/sec 
Rn = (L*V)/v;

Cf = (0.075)/((log(Rn)-2)^2); % Viscous Resistance
displaced_volume = 0.015;
K = 119 * (displaced_volume/(1*B*T) * (B/1))^2;

Cv = Cf + K*Cf;

Rv = Cv*0.5 * rho * V^2 * 0.015;
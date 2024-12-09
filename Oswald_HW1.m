clear;
clc;

%% Problem 1

gamma = 1.4;
M1 = linspace(1,20);
P1 = 1;
T1 = 298.15;

for i = 1:100
M2(i) = sqrt((1+((gamma-1)/2)*M1(i)^2)/(gamma*M1(i)^2-(gamma-1)/2));
P2_over_P1(i) = 1 + 2*gamma/(gamma+1) * (M1(i)^2-1);
T2_over_T1(i) = (1 + 2*gamma/(gamma+1) * (M1(i)^2-1)) * ((2+(gamma-1)*M1(i)^2)/((gamma+1)*M1(i)^2));
rho_ratio(i) = ((gamma+1) * M1(i)^2)/(2+(gamma-1)*M1(i)^2);
total_ratio(i) = M1(i)/M2(i) * ((2 + (gamma-1)*M2(i)^2)/(2 + (gamma-1)*M1(i)^2))^((gamma+1)/(2*(gamma-1)));
end

figure(1)
plot(M1, M2)
title("Problem 1: Mach Number 1 v Mach Number 2")
xlabel("Mach Number 1")
ylabel("Mach Number 2")

figure(2)
plot(M1, P2_over_P1)
title("Problem 1: Mach Number 1 v Pressure Ratio")
xlabel("Mach Number 1")
ylabel("P2/P1")

figure(3)
plot(M1, T2_over_T1)
title("Problem 1: Mach Number 1 v Temperature Ratio")
xlabel("Mach Number 1")
ylabel("T2/T1")

figure(4)
plot(M1, total_ratio)
title("Problem 1: Mach Number 1 v Total Pressure Ratio")
xlabel("Mach Number 1")
ylabel("P02/P01")

figure(5)
plot(M1, rho_ratio)
title("Problem 1: Mach Number 1 v Density Ratio")
xlabel("Mach Number 1")
ylabel("ρ2/ρ1")


%% Problem 2

theta = linspace(0, 120, 121);
M1 = 1;
nu1 = sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(M1^2-1)))-atand(sqrt(M1^2-1));
for i = 1:121
testM = 0.1;
testnu = -999;
nu2(i) = theta(i) + nu1;
while (abs(testnu-nu2(i))) > 0.01
testM = testM + .0001;
testnu = sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)/(gamma+1))*(testM^2-1)))-atand(sqrt(testM^2-1));
end
M2(i) = testM;
T_ratio(i) = (1 + (gamma-1)/2 * M1^2)/(1 + (gamma-1)/2 * M2(i)^2);
P_ratio(i) = T_ratio(i)^(gamma/(gamma-1)); 
end

figure(6)
plot(theta, log10(M2))
title("Problem 2: Turning Angle v log10(Mach Number 2)")
xlabel("Turning Angle(Degrees)")
ylabel("log10(M2)")

figure(7)
plot(theta, log10(P_ratio))
title("Problem 2: Turning Angle v log10(Pressure Ratio)")
xlabel("Turning Angle(Degrees)")
ylabel("log10(P2/P1)")

figure(8)
plot(theta, log10(T_ratio))
title("Problem 2: Turning Angle v log10(Temperature Ratio)")
xlabel("Turning Angle(Degrees)")
ylabel("log10(T2/T1)")

%%Problem 3

Mach1 = linspace(2,50);
theta = 10;

for i = 1:100
testtheta = -999;
testB = 0;

while (abs(testtheta-theta)) > 0.01
testtheta = atand(2*cotd(testB) * ((Mach1(i)^2*sind(testB)^2-1)/(Mach1(i)^2*(gamma+cosd(2*testB))+2)));
testB = testB + .0001;
end
B(i) = testB;
Mn1(i) = Mach1(i) * sind(B(i));
Mn2(i) = sqrt((1+((gamma-1)/2)*Mn1(i)^2)/(gamma*Mn1(i)^2-(gamma-1)/2));
Mach2(i) = Mn2(i) / (sind(B(i) - theta));
P_rat(i) = 1 + 2*gamma/(gamma+1) * (Mn1(i)^2-1);
T_rat(i) = (1 + 2*gamma/(gamma+1) * (Mn1(i)^2-1)) * ((2+(gamma-1)*Mn1(i)^2)/((gamma+1)*Mn1(i)^2));
end


figure(9)
plot(Mach1, Mach2)
title("Problem 3: Mach Number 1 v Mach Number 2")
xlabel("Mach Number 1")
ylabel("Mach Number 2")

figure(10)
plot(Mach1, P_rat)
title("Problem 3: Mach Number 1 v Pressure Ratio")
xlabel("Mach Number 1")
ylabel("Pressure Ratio")

figure(11)
plot(Mach1, T_rat)
title("Problem 3: Mach Number 1 v Temperature Ratio")
xlabel("Mach Number 1")
ylabel("Temperature Ratio")

figure(12)
plot(Mach1, B)
title("Problem 3: Mach Number 1 v Shock Angle")
xlabel("Mach Number 1")
ylabel("Shock Angle(Degrees)")

figure(13)
plot(Mach1, B/theta)
title("Problem 3: Mach Number 1 v Shock/Wedge angle ratio")
xlabel("Mach Number 1")
ylabel("Shock/Wedge Angle Ratio")

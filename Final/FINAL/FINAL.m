%Kathia Coronado
%12/11/20
%RBE 501 Final
%% Problem 1a: solve for dynamical model using lagrange 
clc;
clear;
syms m1 m2 m3 
syms q1 q2 q3 real
syms a b c
syms q_dot_1 q_dot_2 q_dot_3 q_ddot_1 q_ddot_2 q_ddot_3 real
syms g real;

%m2=0;

Slist = [[0;0;1;0;0;0],[0;0;0;0;1;0],[0;0;0;0;0;1]];
thetalist= [q1;q2;q3];
q = [q1; q2; q3;];
q_d = [q_dot_1; q_dot_2; q_dot_3;];
q_dd = [q_ddot_1; q_ddot_2; q_ddot_3;];


 % space
I = eye(3);
T ={};
for i = 1 : length(thetalist)
 w = vector_2_skew(Slist(1:3,i));
 v = Slist(4:6,i);
 theta = thetalist(i);

 R = I + sin(theta)* w + (1-cos(theta))*w^2;

 star = (I * theta+(1-cos(theta))*w+(theta-sin(theta))*w^2)*v;

 T{i} = [R star; 0 0 0 1];

end

M01=[ 0 1 0 0;0 0 1 0;1 0 0 a; 0 0 0 1];
T01=T{1}*M01;

M02=[0 1 0 0;1 0 0 b; 0 0 -1 a;0 0 0 1];
T02=simplify(T{1}*T{2}*M02);
 
M03=[0 1 0 0;1 0 0 b; 0 0 -1 (a-c);0 0 0 1];
T03=simplify(T{1}*T{2}*T{3}*M03);

o1=T01(1:3,4)
o2=T02(1:3,4)
o3=T03(1:3,4)

% velocity kinematics 

Vel_1 = chainRule(o1,q,q_d,q_dd)
Vel_2 = chainRule(o2,q,q_d,q_dd)
Vel_3 = chainRule(o3,q,q_d,q_dd)

%Find lagrangian 
K1 = 0.5 * m1 * Vel_1.' * Vel_1;
K2 = 0.5 * m2 * Vel_2.' * Vel_2;
K3 = 0.5 * m3 * Vel_3.' * Vel_3;
Kin_energy = K1 + K2 + K3;

P1 = m1 * g * T01(3,4);
P2 = m2 * g * T02(3,4);
P3 = m3 * g * T03(3,4);
Pot_energy = P1 + P2 + P3;

L = simplify(Kin_energy - Pot_energy);

%lagrangian equation 
x1 = diff(L, q_dot_1);
x2 = diff(L, q_dot_2);
x3 = diff(L, q_dot_3);

y1 = chainRule(x1, q, q_d, q_dd);
y2 = chainRule(x2, q, q_d, q_dd);
y3 = chainRule(x3, q, q_d, q_dd);

z1 = diff(L, q1);
z2 = diff(L, q2);
z3 = diff(L, q3);

Tau_1 = y1-z1
Tau_2 = y2-z2
Tau_3 = y3-z3

Tau = [Tau_1;Tau_2;Tau_3];

%now get tau in compact form 
M11 = simplify(Tau_1 - subs(Tau_1,q_ddot_1,0)) /q_ddot_1;
M12 = simplify(Tau_1 - subs(Tau_1,q_ddot_2,0)) /q_ddot_2;
M13 = simplify(Tau_1 - subs(Tau_1,q_ddot_3,0)) /q_ddot_3;
M21 = simplify(Tau_2 - subs(Tau_2,q_ddot_1,0)) /q_ddot_1;
M22 = simplify(Tau_2 - subs(Tau_2,q_ddot_2,0)) /q_ddot_2;
M23 = simplify(Tau_2 - subs(Tau_2,q_ddot_3,0)) /q_ddot_3;
M31 = simplify(Tau_3 - subs(Tau_3,q_ddot_1,0)) /q_ddot_1;
M32 = simplify(Tau_3 - subs(Tau_3,q_ddot_2,0)) /q_ddot_2;
M33 = simplify(Tau_3 - subs(Tau_3,q_ddot_3,0)) /q_ddot_3;
M = [[M11,M12,M13]; [M21,M22,M23] ; [M31,M32,M33]];
M = simplify(M)

G1 = subs(Tau_1, [q_dot_1,q_dot_2,q_dot_3,q_ddot_1,q_ddot_2,q_ddot_3], [0,0,0,0,0,0]);
G2 = subs(Tau_2, [q_dot_1,q_dot_2,q_dot_3,q_ddot_1,q_ddot_2,q_ddot_3], [0,0,0,0,0,0]);
G3 = subs(Tau_3, [q_dot_1,q_dot_2,q_dot_3,q_ddot_1,q_ddot_2,q_ddot_3], [0,0,0,0,0,0]);
G = [G1;G2;G3];
G = simplify(G)

C1 = Tau_1 - (M(1,:) * [q_ddot_1,q_ddot_2,q_ddot_3].' + G1);
C2 = Tau_2 - (M(2,:) * [q_ddot_1,q_ddot_2,q_ddot_3].' + G2);
C3 = Tau_3 - (M(3,:) * [q_ddot_1,q_ddot_2,q_ddot_3].' + G3);
C = [C1;C2;C3];
C = simplify(C)

compact_Lagrange_Tau = simplify(expand(M*[q_ddot_1; q_ddot_2; q_ddot_3] + C.*[q_dot_1; q_dot_2; q_dot_3] + G));
clt= subs(compact_Lagrange_Tau, [m2], [0])


%% Problem 1b - derive dynamical model using Newton's approach 

syms I1 I2 I3
syms lc1 lc2 lc3

i = [1; 0; 0];
j = [0; 1; 0];
k = [0; 0; 1];

w0 = [0;0;0];
alpha0 = [0;0;0];
ac0 = [0;0;0];
ae0 = [0;0;0];

R01=[cos(q1) -sin(q1) 0; sin(q1) cos(q1) 0; 0 0 1];
R12=[0 1 0;0 0 1;1 0 0]*[1 0 0;0 1 0 ; 0 0 1];
R02=R01*R12;
R23=[0 0 -1;0 1 0 ; 1 0 0]*[1 0 0;0 1 0 ; 0 0 1];
R03=R02*R23;

z0 = k;
z1 = R02(1:3,3); 
z2 = R03(1:3,3); 
b1 = R01.'*z0;
b2 = R02.'*z1;
b3 = R03.'*z2;
w1 = (R01.'*w0)+(b1*q_dot_1);
w2 = (R12.'*w1)+(b2*q_dot_2);
w3 = (R23.'*w2)+(b3*q_dot_3);
alpha1 = (R01.'*alpha0)+(b1*q_ddot_1)+(cross(w1,b1*q_dot_1));
alpha2 = (R12.'*alpha1)+(b2*q_ddot_2)+(cross(w2,b2*q_dot_2));
alpha3 = (R23.'*alpha2)+(b3*q_ddot_3)+(cross(w3,b3*q_dot_3));

%forward recursion

r1c1 = lc1*k;
r2c1 = (lc1-a)*k;
r12 = a*k;

r2c2 = lc2*j;
r3c2 = (lc2-b)*j;
r23 = b*j;

r3c3 = lc3*-k;
r4c3 = (lc3-(c))*k;


%{
r1c1 = lc1*k;
r2c1 = (lc1-a)*k;
r12 = a*k;

r2c2 = lc2*j;
r3c2 = (lc2-(b))*j;
r23 = (b)*j;

r3c3 = lc3*k;
r4c3 = (lc3-(c+q3))*k;
%}

ac1 = (R01.'*ae0) + (cross(alpha1,r1c1)) + (cross(w1,cross(w1,r1c1)));
ae1 = (R01.'*ae0) + (cross(alpha1,r12)) + (cross(w1,cross(w1,r12)));
ac2 = (R12.'*ae1) + (cross(alpha2,r2c2)) + (cross(w2,cross(w2,r2c2)));
ae2 = (R12.'*ae1) + (cross(alpha2,r23)) + (cross(w2,cross(w2,r23)));
ac3 = (R23.'*ae2) + (cross(alpha3,r3c3)) + (cross(w3,cross(w3,r3c3)));

newton_g1 = -R01.'*g*k;
newton_g2 = -R02.'*g*k;
newton_g3 = -R03.'*g*k;

%backwards recursion

f4 = [0;0;0];
tau4 = [0;0;0];
R_3_4 = eye(3);
f3 = R_3_4*f4 + m3*ac3 - m3*newton_g3;
tau3 = R_3_4*tau4 - cross(f3,r3c3) + cross(R_3_4*f4,r4c3) + I3*alpha3 + cross(w3,I3*w3);
f2 = R23*f3 + m2*ac2 - m2*newton_g2;
tau2 = R23*tau3 - cross(f2,r2c2) + cross(R23*f3,r3c2) + I2*alpha2 + cross(w2,I2*w2);
f1 = R12*f2 + m1*ac1 - m1*newton_g1;
tau1 = R12*tau2 - cross(f1,r1c1) + cross(R12*f2,r2c1) + I1*alpha1 + cross(w1,I1*w1);
Newton_Tau = simplify(expand([tau1(3); tau2(3); tau3(3)]));

M11 = simplify(tau1(3) - subs(tau1(3),q_ddot_1,0)) /q_ddot_1;
M12 = simplify(tau1(3) - subs(tau1(3),q_ddot_2,0)) /q_ddot_2;
M13 = simplify(tau1(3) - subs(tau1(3),q_ddot_3,0)) /q_ddot_3;
M21 = simplify(tau2(3) - subs(tau2(3),q_ddot_1,0)) /q_ddot_1;
M22 = simplify(tau2(3) - subs(tau2(3),q_ddot_2,0)) /q_ddot_2;
M23 = simplify(tau2(3) - subs(tau2(3),q_ddot_3,0)) /q_ddot_3;
M31 = simplify(tau3(3) - subs(tau3(3),q_ddot_1,0)) /q_ddot_1;
M32 = simplify(tau3(3) - subs(tau3(3),q_ddot_2,0)) /q_ddot_2;
M33 = simplify(tau3(3) - subs(tau3(3),q_ddot_3,0)) /q_ddot_3;
Mn = [[M11,M12,M13]
 [M21,M22,M23]
 [M31,M32,M33]];

G1 = subs(tau1(3), [q_dot_1,q_dot_2,q_dot_3,q_ddot_1,q_ddot_2,q_ddot_3], [0,0,0,0,0,0]);
G2 = subs(tau2(3), [q_dot_1,q_dot_2,q_dot_3,q_ddot_1,q_ddot_2,q_ddot_3], [0,0,0,0,0,0]);
G3 = subs(tau3(3), [q_dot_1,q_dot_2,q_dot_3,q_ddot_1,q_ddot_2,q_ddot_3], [0,0,0,0,0,0]);
Gn = [G1;G2;G3];
Gn = simplify(Gn);


C1 = tau1(3) - (Mn(1,:) * [q_ddot_1,q_ddot_2,q_ddot_3].' + G1);
C2 = tau2(3) - (Mn(2,:) * [q_ddot_1,q_ddot_2,q_ddot_3].' + G2);
C3 = tau3(3) - (Mn(3,:) * [q_ddot_1,q_ddot_2,q_ddot_3].' + G3);
Cn = [C1;C2;C3];
Cn = simplify(Cn);

compact_Newton_Tau = simplify(expand(Mn*[q_ddot_1; q_ddot_2; q_ddot_3] + Cn.*[q_dot_1; q_dot_2; q_dot_3] + Gn))

cnt= subs(compact_Newton_Tau,[I1,I2,I3,lc1,lc2,lc3,m2],[0,0,0,a,0,c+q3,0])

%% Problem 1c 

a_=0.3;
b_=0.1;
c_=0.1;
m1_=0.5;
m3_=0.5;
m2_=0;
g_=9.8;

Lan= subs(compact_Lagrange_Tau,[a,b,c,m1,m2,m3,g],[a_,b_,c_,m1_,m2_,m3_,g_])
Newt= subs(compact_Newton_Tau,[a,b,c,m1,m2,m3,g],[a_,b_,c_,m1_,m2_,m3_,g_])


%% Problem 1d 
clc
clear all
close all

syms q1 q2 q3
syms t

a = 0.3;
b = 0.1;
c = 0.1;
p1 = [0; 0; a];
p2 = [-sin(q1)*(b + q2); cos(q1)*(b + q2); a];
pe = [-sin(q1)*(b + q2); cos(q1)*(b + q2); a-c+q3];
a0 = [0; 200; 150];
a1 = [0; 0; 0];
a2 = [6; -6; 1.5];
a3 = [-0.4; 0.4; -0.1];
x = 0.01:0.01:10;
p = a0 + (a1*t) + (a2*(t^2)) + (a3*(t^3))
p_dot = a1 + (2 * a2 * (t)) + (3 * a3 * (t^2))
p_dot_dot = (2 * a2) + (6 * a3 * (t))

straight_task_set = subs(p, t, x);
figure()
plot3(straight_task_set(1,:), straight_task_set(2,:), straight_task_set(3,:));

title ('position')
xlabel('position (mm)')
ylabel('position (mm)')
zlabel('position (mm)')
axis equal

straight_task_set = subs(p_dot, t, x);
figure()
plot3(straight_task_set(1,:), straight_task_set(2,:), straight_task_set(3,:));

title ('velocity')
xlabel('velocity (mm/s)')
ylabel('velocity (mm/s)')
zlabel('velocity (mm/s)')
axis equal

straight_task_set = subs(p_dot_dot, t, x);
figure()
plot3(straight_task_set(1,:), straight_task_set(2,:), straight_task_set(3,:));

title ('accn')
xlabel('accn (mm/s^2)')
ylabel('accn (mm/s^2)')
zlabel('accn (mm/s^2)')
axis equal

eqn_straight_line_for_q = p == pe;

[q1_sol q2_sol q3_sol] = solve(eqn_straight_line_for_q, [q1 q2 q3])
double(subs([q1_sol q2_sol q3_sol], t, 0.01))
q1_sol
q2_sol
q3_sol

straight_q1_set = subs(q1_sol(2), t, x);
straight_q2_set = subs(q3_sol(2), t, x);
straight_q3_set = subs(q3_sol(2), t, x);
figure()
plot(x, straight_q1_set);
hold on
plot(x, straight_q2_set);
hold on
plot(x, straight_q3_set);
title('joint-position plot');
legend('joint1', 'joint2', 'joint3');
figure()
plot(x, straight_q1_set);
title('1st joint time-position plot');
xlabel('time (s)')
ylabel('joint-position')
figure()
plot(x, straight_q2_set);
title('2nd joint time-position plot');
xlabel('time (s)')
ylabel('joint position')
figure()
plot(x, straight_q3_set);
title('3rd joint time-position plot');
xlabel('time (s)')
ylabel('joint position')

q1_dot = diff(q1_sol(2), t);
q2_dot = diff(q2_sol(2), t);
q3_dot = diff(q3_sol(2), t);
straight_q1_dot_set = subs(q1_dot, t, x);
straight_q2_dot_set = subs(q2_dot, t, x);
straight_q3_dot_set = subs(q3_dot, t, x);
figure()
plot(x, straight_q1_dot_set);
hold on
plot(x, straight_q2_dot_set);
hold on
plot(x, straight_q3_dot_set);
title('joint-velocity plot');
legend('joint1', 'joint2', 'joint3');
figure()
plot(x, straight_q1_dot_set);
title('1st joint time-velocity plot');
xlabel('time (s)')
ylabel('joint velocity')
figure()
plot(x, straight_q2_dot_set);
title('2nd joint time-velocity plot');
xlabel('time (s)')
ylabel('joint velocity')
figure()
plot(x, straight_q3_dot_set);
title('3rd joint time-velocity plot');
xlabel('time (s)')
ylabel('joint velocity')

q1_ddot = diff(diff(q1_sol(2), t), t);
q2_ddot = diff(diff(q2_sol(2), t), t);
q3_ddot = diff(diff(q3_sol(2), t), t);
straight_q1_ddot_set = subs(q1_ddot, t, x);
straight_q2_ddot_set = subs(q2_ddot, t, x);
straight_q3_ddot_set = subs(q3_ddot, t, x);
figure()
plot(x, straight_q1_ddot_set);
hold on
plot(x, straight_q2_ddot_set);
hold on
plot(x, straight_q3_ddot_set);
title('joint acceleration plot');
legend('joint1', 'joint2', 'joint3');
figure()
plot(x, straight_q1_ddot_set);
title('1st joint time-accn plot');
xlabel('time (s)')
ylabel('joint acc')
figure()
plot(x, straight_q2_ddot_set);
title('2nd joint time-accn plot');
xlabel('time (s)')
ylabel('joint acc')
figure()
plot(x, straight_q3_ddot_set);
title('3rd joint time-accn plot');
xlabel('time (s)')
ylabel('joint acc')


xyz = subs(pe, [q1, q2, q3], [q1_sol(2), q2_sol(2), q3_sol(2)]);
xyz_set = subs(xyz, t, x);
figure()
plot(x, xyz_set(3,:));
hold on
plot(x, straight_q3_set);
title('Verify')
axis equal
legend('Solution trajectory','Desired Trajectory')

%%

function answer = chainRule(func,q,qd,qdd)
    for i = 1:length(func)
        for j = 1:length(q)                
            partials(:,j) = diff(func(i),q(j)) * qd(j) + diff(func(i),qd(j)) * qdd(j);
        end             
        temp(:,i) = simplify(partials);
    end       
    answer = sum(temp).';
end

function X = vector_2_skew(x)

X=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];

end

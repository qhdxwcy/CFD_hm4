%% 第三题作图
x = 0:0.01:3.2;
y1 = sin(x);
y2 = 2*sin(x)-1/2*sin(2*x);
y3 = 4/3*sin(x)-1/6*sin(2*x);
plot(x,y1);
hold on;
grid on;
plot(x,y2);
plot(x,y3);
axis equal;
xlim([0,3.2]);
xlabel('k');
ylabel('Re(k'')');
legend('一阶迎风','二阶迎风','三阶偏迎风','Location','best');
hold off;

x = 0:0.01:3.5;
z1 = cos(x)-1;
z2 = 2*cos(x)-1/2*cos(2*x)-3/2;
z3 = 2/3*cos(x)-1/6*cos(2*x)-1/2;
plot(x,z1);
hold on;
grid on;
plot(x,z2);
plot(x,z3);
xlim([0,3.5]);
ylim([-3,0]);
axis equal;
xlabel('k');
ylabel('Im(k'')');
legend('一阶迎风','二阶迎风','三阶偏迎风','Location','best');
hold off;

%% 计算分辨能力
%一阶迎风
x = 0:0.01:pi;
y1 = abs((sin(x)-x)./x)-0.005;
y2 = abs((cos(x)-1-x)./x)-0.005;
plot(x,min(y1,y2))

f = @(k) abs((sin(k)-k)./ k) - 0.005;
g = @(k) abs((cos(k)-1-k) ./ k) - 0.005;
k_guess = 0.2;
options = optimset('Display', 'off');
k_solution1 = fsolve(f, k_guess, options)
k_solution2 = fsolve(g, k_guess, options)

%二阶迎风
f = @(k) abs((2*sin(k)-1/2*sin(2*k)-k)./ k) - 0.005;
g = @(k) abs((2*cos(k)-1/2*cos(2*k)-3/2-k) ./ k) - 0.005;
k_guess = 0.2;
options = optimset('Display', 'off');
k_solution3 = fsolve(f, k_guess, options)
k_solution4 = fsolve(g, k_guess, options)

%三阶偏迎风
f = @(k) abs((4/3*sin(k)-1/6*sin(2*k)-k)./ k) - 0.005;
g = @(k) abs((2/3*cos(k)-1/6*cos(2*k)-1/2-k) ./ k) - 0.005;
k_guess = 0.2;
options = optimset('Display', 'off');
k_solution5 = fsolve(f, k_guess, options)
k_solution6 = fsolve(g, k_guess, options)


%% 初始条件
Mx=100;%空间网格数
Nt=2000;%时间网格数
N=Mx*Nt; %总网格数
delta_x = 1/Mx;
delta_t = 10/Nt;
c = 1*delta_t/delta_x

U=zeros(1,Mx+1);
%先赋初值和边界值
for x = 1:Mx+1
    U(x)=InitialConditions(-0.5+(x-1)*delta_x);
end
U(Mx+1)=U(1);%边界周期条件


%% Lax-Wendroff格式
M = zeros(1,Mx+1);
UU = zeros(3,Mx+1);
T = [Nt/100 Nt/10 Nt];
for i =1:3
    for t=1:T(i)
        for x=2:Mx
            M(x) = U(x)+0.5*c*(U(x-1)-U(x+1))+0.5*c^2*(U(x-1)-2*U(x)+U(x+1));
        end
        M(1) = U(1)+0.5*c*(U(Mx)-U(2))+0.5*c^2*(U(Mx)-2*U(1)+U(2));
        M(Mx+1)=U(1);
        U = M;
    end
    UU(i,:) = U;
end

%% Warming-Beam格式
M = zeros(1,Mx+1);
UU = zeros(3,Mx+1);
T = [Nt/100 Nt/10 Nt];
for i =1:3
    for t=1:T(i)
        for x=3:Mx
            M(x) = U(x)+0.5*c*(4*U(x-1)-3*U(x)-U(x-2))+0.5*c^2*(U(x-2)-2*U(x-1)+U(x));
        end
        M(1) = U(1)+0.5*c*(4*U(Mx)-3*U(1)-U(Mx-1))+0.5*c^2*(U(Mx-1)-2*U(Mx)+U(1));
        M(2) = U(2)+0.5*c*(4*U(1)-3*U(2)-U(Mx))+0.5*c^2*(U(Mx)-2*U(1)+U(2));
        M(Mx+1)=U(1);
        U = M;
    end
    UU(i,:) = U;
end

%% 新格式
M = zeros(1,Mx+1);
UU = zeros(3,Mx+1);
T = [Nt/100 Nt/10 Nt];
for i =1:3
    for t=1:T(i)
        for x=3:Mx
            M(x) = U(x)-1/6*c*((2-c)*U(x+1)-(6+3*c)*U(x-1)+(3+3*c)*U(x)+(1+c)*U(x-2))+1/6*c^2*((1+c)*U(x-2)-3*c*U(x-1)+(3*c-3)*U(x)+(2-c)*U(x+1));
        end
        M(1) = U(1)-1/6*c*((2-c)*U(2)-(6+3*c)*U(Mx)+(3+3*c)*U(1)+(1+c)*U(Mx-1))+1/6*c^2*((1+c)*U(Mx-1)-3*c*U(Mx)+(3*c-3)*U(1)+(2-c)*U(2));
        M(2) = U(2)-1/6*c*((2-c)*U(3)-(6+3*c)*U(1)+(3+3*c)*U(2)+(1+c)*U(Mx))+1/6*c^2*((1+c)*U(Mx)-3*c*U(x-1)+(3*c-3)*U(2)+(2-c)*U(3));
        M(Mx+1)=U(1);
        U = M;
    end
    UU(i,:) = U;
end


%% 作图
pA =  UU(1,:);    
pB =  UU(2,:);
pC =  UU(3,:);
XX = -0.5:delta_x:0.5;
plot(XX,pA)
hold on;
plot(XX,pB)
hold on;
plot(XX,pC)
xlabel('x');
ylabel('u');
legend('t=0.1','t=1.0','t=10.0','Location','best');
title('Warming-Beam格式数值解');


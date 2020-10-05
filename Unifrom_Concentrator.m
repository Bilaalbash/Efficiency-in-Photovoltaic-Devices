%% ---- Project Information ----%
% Student Name: Bilaal Bashir
% Student Number: 201199049
% Supervisor: Dr Zoran Ikonic
% Project Name: Efficient Photovoltaic Devices for Renewable Energy Sources Using Solar Concentrators

%% ---- Clear Command Window ----%    

clear all
clc 
close all 
format long

%% ---- Setting Paremeters ----%    

y_values = [-10 -40 -70];

x(1) = 10; %x0
N = 3.882 + 0.019i;
N1 = 2.757 + 3.792i;

RGI = 0.613; %reflection
TGI = 0.6515; %transmittance
C = 1000; %concentration ratio

DX = 0.01; %step size

%% ---- Setting y(1) for plotting ----%    
for ii=1:1:3

    if ii == 1
        y(1) = y_values(1);
    elseif ii == 2
        y(1) = y_values(2);
    elseif ii == 3
        y(1) = y_values(3);
        
    end
    
    g(1) = (sqrt((y(1))^2+(x(1))^2)+y(1))/(x(1)); %g0

%% for loop

    for i=1:1:10001; % 1 to 100 in steps of 1
 %% ---- Calculating x(i) ----%    


        x_05 = x(i) + ((DX)/(2)); % x(i+0.5)

        x(i+1) = x_05 + ((DX)/(2)); % x(i+1)
    
%% ---- Calculating R(g(i)) ----%    

        a(i+1) = sqrt((N1^2-1)*(g(i)^2) + N1^2) - 1;

        b(i+1) = sqrt((N1^2-1)*(g(i)^2) + N1^2) + 1;

        c(i+1) = sqrt((N1^2-1)*(g(i)^2) + N1^2) - g(i)^2;

        d(i+1) = sqrt((N1^2-1)*(g(i)^2) + N1^2) + g(i)^2;


        R(i+1) = ((0.5)*(((abs(a(i)/b(i)))^2) * (1 + (abs(c(i)/d(i)))^2)));
    
%% ---- Calculating T(g(i)) ----%    

        l(i+1) = (sqrt((N^2*(1+g(i)^2)^2) - (1-g(i)^2)^2) - (2*g(i))); 

        m(i+1) = (sqrt((N^2*(1+g(i)^2)^2) - (1-g(i)^2)^2) + (2*g(i)));

        n(i+1) = (2*g(i))*(sqrt((N^2*(1+g(i)^2)^2) - (1-g(i)^2)^2)) - (1-(g(i))^2)^2;

        o(i+1) = (2*g(i))*(sqrt((N^2*(1+g(i)^2)^2) - (1-g(i)^2)^2)) + (1-(g(i))^2)^2;

    % Transmitivity
        T(i+1) = 1 - (0.5)*((abs(l(i)/m(i)))^2) * (1 + (abs(n(i)/o(i)))^2);
    
%% ---- Implementing Runge Kutta ----%      
        
        Function = (g(i)*(g(i)^2+1)*(RGI)*(TGI)*(C)-(2*g(i)^2))/(x(i)*(g(i)^2+1)*(RGI)*(TGI)*(C)); 

        g_05 = g(i) + ((DX/2) * Function);
    
        Function_05 = (g_05*(g_05^2+1)*(RGI)*(TGI)*(C)-(2*g_05^2))/(x_05*(g_05^2+1)*(RGI)*(TGI)*(C));
   
        g(i + 1) = g(i) + (DX * Function_05);
    
        y(i+1) = y(1) + (g(i)* x(i)/2);
    
    end
    
    if ii == 1
        R_10 = R;
        T_10 = T;
        y_10 = y;
    elseif ii == 2
        R_40 = R;
        T_40 = T;
        y_40 = y;
    elseif ii == 3
        R_70 = R;
        T_70 = T;
        y_70 = y;
        
    end

end

%% ---- Figures ----%    

%Reflection against x (cm)
%figure
subplot(2,2,2)
plot(x,R_10,'r','LineWidth',1.1)
hold on
plot(x,R_40,'b','LineWidth',1.1)
plot(x,R_70,'g','LineWidth',1.1)
hold off
set(gca,'XLim',[10 110])
set(gca,'YLim',[0.56 0.62])
xlabel('x (cm)')
ylabel('R')
legend({'y0 = -10 cm','y0 = -40 cm', 'y0 = -70 cm'}, 'Location', 'southwest')
grid


%Transmittivity against x (cm)
%figure
subplot(2,2,3)
plot(x,T_10,'r','LineWidth',1.1)
hold on
plot(x,T_40,'b','LineWidth',1.1)
plot(x,T_70,'g','LineWidth',1.1)
hold off
set(gca,'XLim',[10 110])
set(gca,'YLim',[0.52 0.66])
xlabel('x (cm)')
ylabel('T')
legend({'y0 = -10 cm','y0 = -40 cm', 'y0 = -70 cm'}, 'Location', 'southeast')
grid;

%Reflection Profile y(x)
%figure
subplot(2,2,1)
plot(x,y_10,'r','LineWidth',1.1)
hold on
plot(x,y_40,'b','LineWidth',1.1)
plot(x,y_70,'g','LineWidth',1.1)
hold off
set(gca,'XLim',[10 110])
set(gca,'YLim',[-100 250])
xlabel('x (cm)')
ylabel('y(x) (cm)')
legend({'y0 = -10 cm','y0 = -40 cm', 'y0 = -70 cm'}, 'Location', 'northwest')
grid;


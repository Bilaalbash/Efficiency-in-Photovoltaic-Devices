%% ---- Project Information ----%
% Student Name: Bilaal Bashir
% Student Number: 201199049
% Supervisor: Dr Zoran Ikonic
% Project Name: Efficient Photovoltaic Devices for Renewable Energy Sources Using Solar Concentrators

%% --- Clear Command Window --- %%
clear all
clc 
close all 
format long

%% --- Setting Paremeters ---%    
R = 1;
r_values = [1 0.8 0.5 0.2];
Sum = 0;
dash_sum = 0;
a_dash_temp = 0;
k = 1;
m = 1;
A_final = 0;
A_final_save = 0;
circular_final = 0;
 
for ik = 1:4
     
    if ik == 1 
          r(1) = r_values(1);
    elseif ik == 2
          r(1) = r_values(2);
    elseif ik == 3
          r(1) = r_values(3);
    elseif ik == 4
          r(1) = r_values(4);
    end
     
    for i = 1:44
     
        theta(i) = 45 + i; % Increments Theta from 46 to 90
        i_max(i) = 45 / (90 - theta(i)); % Max number of reflections
        i_int(i) = floor(i_max(i)); % Stores the integer value of i_max
        f(i) = i_max(i) - i_int(i); % Stores the float value of i_max

    % Calculating l_max depending on float value of i_max
    
        if f(i) == 0 
            l_max(i) = i_int(i) - 1; %l_max is the max number of layers
    
        elseif f(i) > 0
            l_max(i) = i_int(i);
        
        end

    % Calculating W_max depending on float value of i_max

        if f(i) == 0
            W_max(i) = 2*R; %W_max is the width at the maximum layer (l_max)
    
        elseif f(i) > 0
            W_max(i) = 2*f(i)*R;
        end

        P_max(i) = 90 + ((l_max(i) + 1) - 1)*(2*theta(i) - 180); % P_lmax+1
        a_max(i) = ((W_max(i)) * tand(P_max(i))) / (tand(theta(i)) - tand(P_max(i))); % a when l = l_max
        a_sum(i) = a_max(i); % the first value of a is saved into the sum of a
    
        A(1) = a_sum(i);
        W(1) = W_max(i);
    
        for j=1:100
    
            if l_max(i) > 1 % if l_max is greater than one 
            
                l(j) = l_max(i) - j; % decrease l_max by 1 each loop
                W_run = 2*(R + a_sum(i));
            
            elseif l_max(i) == 1  % if l_max is equal to one 
                break % break the loop
            end
    
            P_run = 90 + ((l(j) + 1) - 1)*(2*theta(i) - 180); % Calculate new value of P
            a_run = ((W_run) * tand(P_run)) / (tand(theta(i)) - tand(P_run)); % Calculate new value of a
            a_sum(i) = a_run + a_sum(i); % add new value of a to the sum of a
        
            A(j+1) = a_run; % saves all values of a into an array
            W(j+1) = W_run;% saves all values of W into an array

        %% --- Calculation a' values ---%%

            if l(j)<=1 % if l is less than or equal to one break the loop
                break
            end

        end

        w = flip(W); %flips the array so the biggest value of W is first for simplicity in future calculations
        a = flip(A); %flips the array so the biggest value of a is first for simplicity in future calculations
        Sum = 0;
        dash_sum = 0;
        a_dash_temp = 0;
        A_final = 0;
        circular_final = 0;

        for t = 1:1:l_max(i)
       
            Sum = Sum + a(t);
            a_dash_temp = ((Sum/w(t)) * (2*R)) - dash_sum;
       
            if t == l_max(i)
                a_dash_temp = ((Sum) - dash_sum);
            end
            
            
            second = (R + dash_sum)^2;
       
            dash_sum_1 = dash_sum;
            dash_sum = dash_sum + a_dash_temp;
            first = (R + dash_sum)^2;
            dash_sum_save(k) = a_dash_temp;
            
            circular_final = circular_final + (r^t * (first - second));
       
            A_final = A_final + (a_dash_temp * r^t);
            A_final_save(t) = (a_dash_temp * r^t);
       

        end
   
   %c_plus(i) = 1 + 2(*r^l*a_dash_temp);
        if ik == 1
            c_plus(i) = 1 + 2*(A_final/R);
            c_circ(i) = 1 + (circular_final / R^2);
        elseif ik == 2
            c_plus_2(i) = 1 + 2*(A_final/R);
            c_circ_2(i) = 1 + (circular_final / R^2);
        elseif ik == 3
            c_plus_3(i) = 1 + 2*(A_final/R);
            c_circ_3(i) = 1 + (circular_final / R^2);
        elseif ik == 4
            c_plus_4(i) = 1 + 2*(A_final/R);
            c_circ_4(i) = 1 + (circular_final / R^2);
        end
        
        H(i) = sum(A) * tand(theta(i)); % calculates the height of the cone

end
        Sum = 0;
        dash_sum = 0;
        a_dash_temp = 0;
        A_final = 0;
        circular_final = 0;
        clear A;
        clear W;
        clear i_int;
        clear f;

 end
 
 %% --- Plot the Figures --- %%%

subplot(3,2,1)
yyaxis left
plot(theta, i_max)
set(gca,'YLim',[1 100])
set(gca, 'YScale', 'log')
ylabel('i max')
yyaxis right
plot(theta, H)
xlabel('Cone Angle (Degree)')
ylabel('Height of Cone (m)')
set(gca,'XLim',[45 90])
set(gca, 'YScale', 'log')
hold on

subplot(3,2,2)
yyaxis right
plot(theta, H)
xlabel('Cone Angle (Degree)')
ylabel('Height of Cone (m)')
set(gca,'XLim',[45 90])
set(gca,'YLim',[0.001 1250])
set(gca, 'YScale', 'log')
yyaxis left
set(gca,'YLim',[1 7000])
set(gca, 'YScale', 'log')
ylabel('Concentration Ratio')
hold on
plot(theta, c_plus,'b')
plot(theta, c_plus_2,'b') 
plot(theta, c_circ,'b')
plot(theta, c_circ_2,'b')
legend({'Plus (r = 1.0)','Plus (r = 0.8)','Circular (r = 1.0)','Circular (r = 0.8)'}, 'Location', 'northwest')
 
subplot(3,2,3)
xlabel('Cone Angle (Degree)')
set(gca,'YLim',[1 7000])
set(gca, 'YScale', 'log')
ylabel('Concentration Ratio')
hold on
plot(theta, c_plus,'b','LineWidth',1)
plot(theta, c_circ,'r','LineWidth',1)
legend({'Plus (r = 1.0)','Circular (r = 1.0)','Plus (r = 0.8)','Circular (r = 0.8)'}, 'Location', 'northwest')

subplot(3,2,4)
xlabel('Cone Angle (Degree)')
set(gca,'YLim',[1 7000])
set(gca, 'YScale', 'log')
ylabel('Concentration Ratio')
hold on
plot(theta, c_plus_2,'b','LineWidth',1)
plot(theta, c_circ_2,'r','LineWidth',1)
legend({'Plus (r = 0.8)','Circular (r = 0.8)'}, 'Location', 'northwest')

subplot(3,2,5)
xlabel('Cone Angle (Degree)')
set(gca,'YLim',[1 7000])
set(gca, 'YScale', 'log')
ylabel('Concentration Ratio')
hold on
plot(theta, c_plus_3,'b','LineWidth',1)
plot(theta, c_circ_3,'r','LineWidth',1)
legend({'Plus (r = 0.5)','Circular (r = 0.5)'}, 'Location', 'northwest')
 
subplot(3,2,6)
xlabel('Cone Angle (Degree)')
set(gca,'YLim',[1 7000])
set(gca, 'YScale', 'log')
ylabel('Concentration Ratio')
hold on
plot(theta, c_plus_4,'b','LineWidth',1)
plot(theta, c_circ_4,'r','LineWidth',1)
legend({'Plus (r = 0.2)','Circular (r = 0.2)'}, 'Location', 'northwest')
 
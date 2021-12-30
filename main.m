%fitting GNSS deformation time series data for afterslip (exponent + Logarithm) using simmulated annealing
%algorithm (SA)
%Mohammad Rheza Zamani
clear
clc;
%Import data
Us_data = importdata('LEWK-east.txt');
t = importdata('LEWK-date.txt');
%Parameter modelling
nitr = 500; %Jumlah iterasi 
T = 5;
dec = 0.01;
%Search space
A_min = -10;
A_max = 1;
B_min = -200;
B_max = 1;
C_min = 1;
C_max = 10;
D_min = 1;
D_max = 10;
V0_min = 1;
V0_max = 200;
Ta_min = 1;
Ta_max = 50;
Tb_min = 1;
Tb_max = 50;
Tc_min = 1;
Tc_max = 50;
%Initial model
%model1(1,1) = A_min + rand*(A_max-A_min);
%model1(1,2) = B_min + rand*(B_max-B_min);
%model1(1,3) = C_min + rand*(C_max-C_min);
%model1(1,4) = D_min + rand*(D_max-D_min);
%model1(1,5) = V0_min + rand*(V0_max-V0_min);
%model1(1,6) = Ta_min + rand*(Ta_max-Ta_min);
%model1(1,7) = Tb_min + rand*(Tb_max-Tb_min);
%model1(1,8) = Tc_min + rand*(Tc_max-Tc_min);
model1(1,1) = 20;
model1(1,2) = 20;
model1(1,3) = 20;
model1(1,4) = 20;
model1(1,5) = 20;
model1(1,6) = 20;
model1(1,7) = 20;
model1(1,8) = 20;
Us_cal1 = cal_slip(model1(1),model1(2),model1(3),model1(4),model1(5),model1(6),model1(7),model1(8),t);
E1 = fun_obj(Us_data,Us_cal1);
v = VideoWriter('Afterslip modelling Easting.avi');
open(v);
for itr = 1 : nitr
    model2(1,1) = A_min + rand*(A_max-A_min);
    model2(1,2) = B_min + rand*(B_max-B_min);
    model2(1,3) = C_min + rand*(C_max-C_min);
    model2(1,4) = D_min + rand*(D_max-D_min);
    model2(1,5) = V0_min + rand*(V0_max-V0_min);
    model2(1,6) = Ta_min + rand*(Ta_max-Ta_min);
    model2(1,7) = Tb_min + rand*(Tb_max-Tb_min);
    model2(1,8) = Tc_min + rand*(Tc_max-Tc_min);
    Us_cal2 = cal_slip(model2(1),model2(2),model2(3),model2(4),model2(5),model2(6),model2(7),model2(8),t);
    E2 = fun_obj(Us_data,Us_cal2);
    delta_E = E2 - E1;
    if delta_E<0
        model1 = model2;
        E1 = E2;
    else
        P = exp(-delta_E/T); 
        if P>= rand
           model1 = model2;
           E1 = E2;
        end
    end
    Us_new = cal_slip(model1(1),model1(2),model1(3),model1(4),model1(5),model1(6),model1(7),model1(8),t);
   Egen(itr) = E1;
   T = T*(1-dec);
   Temperature(itr) = T;
% Plotting 
figure(1)
scatter(t,Us_data,'b.')
hold on
plot(t,Us_new,'r')
xlabel('Day','FontSize',10,'FontWeight','Bold')
ylabel('Displacement Us (cm)','FontSize',10,'FontWeight','Bold')
title(['Exponential and Logarithm models of Easting LEWK|| Hasil Inversi ke-.',num2str(itr)],'FontWeight','bold')
subtitle(['A = ',num2str(model1(1)),' ; B = ',num2str(model1(2)),' ; C = ',num2str(model1(3)),' ; D = ',num2str(model1(4)),' ; V0 = ',num2str(model1(5)),' ; Ta = ',num2str(model1(6)),' ; Tb = ',num2str(model1(7)),' ; Tc = ',num2str(model1(8)), '|| ERMS = ',num2str(Egen(itr))],'FontWeight','bold')
legend('GPS data','models')
grid on
set(gcf, 'Position', get(0, 'Screensize'));
hold off
frame = getframe(gcf);
writeVideo(v,frame);
end
close(v);
%plot misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on

%Plot Temperature
figure(3)
plot(1:nitr,Temperature,'b','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('Temperature','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Penurunan Temperature ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on

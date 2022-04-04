%fitting GNSS deformation time series data for afterslip (exponent +
%Logarithm) using modified Symbiotic Organism Search
%Mohammad Rheza Zamani
clear
clc;
%Import data
Us_data = importdata('LEWK-east.txt');
t = importdata('LEWK-date.txt');
%Parameter modelling
nitr = 500; %Jumlah iterasi 
npop = 100;
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
for i = 1 : npop
    model(i,1) = A_min + rand*(A_max-A_min);
    model(i,2) = B_min + rand*(B_max-B_min);
    model(i,3) = C_min + rand*(C_max-C_min);
    model(i,4) = D_min + rand*(D_max-D_min);
    model(i,5) = V0_min + rand*(V0_max-V0_min);
    model(i,6) = Ta_min + rand*(Ta_max-Ta_min);
    model(i,7) = Tb_min + rand*(Tb_max-Tb_min);
    model(i,8) = Tc_min + rand*(Tc_max-Tc_min);
    Us_cal(i,:) = cal_slip(model(i,1),model(i,2),model(i,3),model(i,4),model(i,5),model(i,6),model(i,7),model(i,8),t);
    E(i) = fun_obj(Us_data,Us_cal(i,:));
end
v = VideoWriter('Afterslip modelling Easting SOS.avi');
open(v);
%Proses inversi
for itr = 1 : nitr
    for i = 1 : npop
        idx = find(E ==min(E));
        model_best = model(idx(1),:);
        %Mutualisme
        j = randi(npop,1);
        k = randi(npop,1);
        if j==i || k==i
            j = randi(npop,1);
            k = randi(npop,1);
        end
        model_mut = [model(i,:);model(j,:)];
        mv_m =(model(i,:)+model(j,:))/2;
        bf = 1;
        for l = 1 : 2
            mod_mut(l,:) = model_mut(l,:) + rand*(model(k)-mv_m*bf);
            if mod_mut(l,1)<A_min 
                mod_mut(l,1) = A_min;
            end
            if mod_mut(l,2)< B_min
                mod_mut(l,2) = B_min;
            end
            if mod_mut(l,3)< C_min
                mod_mut(l,3) = C_min;
            end
            if mod_mut(l,4)<D_min
                mod_mut(l,4) = D_min;
            end
            if mod_mut(l,5)<V0_min
                mod_mut(l,5) = V0_min;
            end
            if mod_mut(l,6)<Ta_min
                mod_mut(l,6) = Ta_min;
            end
            if mod_mut(l,7) < Tb_min
                mod_mut(l,7) = Tb_min;
            end
            if mod_mut(l,8) < Tc_min
                mod_mut(l,8) = Tc_min;
            end

            if mod_mut(l,1)>A_max 
                mod_mut(l,1) = A_max;
            end
            if mod_mut(l,2)> B_max
                mod_mut(l,2) = B_max;
            end
            if mod_mut(l,3)> C_max
                mod_mut(l,3) = C_max;
            end
            if mod_mut(l,4)>D_max
                mod_mut(l,4) = D_max;
            end
            if mod_mut(l,5)>V0_max
                mod_mut(l,5) = V0_max;
            end
            if mod_mut(l,6)>Ta_max
                mod_mut(l,6) = Ta_max;
            end
            if mod_mut(l,7) > Tb_max
                mod_mut(l,7) = Tb_max;
            end
            if mod_mut(l,8) > Tc_max
                mod_mut(l,8) = Tc_max;
            end
        end
        %Hitung model untuk prosedur mutualisme
        for l = 1 : 2
            [cal_slip_mut] = cal_slip(mod_mut(l,1),mod_mut(l,2),mod_mut(l,3),mod_mut(l,4),mod_mut(l,5),mod_mut(l,6),mod_mut(l,7),mod_mut(l,8),t);
            err_mut(l) = fun_obj(Us_data,cal_slip_mut);
            %Update model jika  nilai misfit lebih baik proses mutualisme
            if l == 1
                if err_mut(l)<E(i)
                    model(i,:) = mod_mut(l,:);
                    E(i) = err_mut(l);
                    Us_cal(i,:) = cal_slip_mut;
                end
            else
                if err_mut(l)<E(j)
                    model(j,:) = mod_mut(l,:);
                    E(j) = err_mut(l);
                    Us_cal(i,:) = cal_slip_mut;
                end
            end
        end
        %Komensalisme
        j = randi(npop,1);
        if j == i
            j = randi(npop,1);
        end
        mod_com = model(i) +(0.4+0.9*rand)*(model_best-model(j));
            if mod_com(1)<A_min 
                mod_com(1) = A_min;
            end
            if mod_com(2)< B_min
                mod_com(2) = B_min;
            end
            if mod_com(3)< C_min
                mod_com(3) = C_min;
            end
            if mod_com(4)<D_min
                mod_com(4) = D_min;
            end
            if mod_com(5)<V0_min
                mod_com(5) = V0_min;
            end
            if mod_com(6)<Ta_min
                mod_com(6) = Ta_min;
            end
            if mod_com(7) < Tb_min
                mod_com(7) = Tb_min;
            end
            if mod_com(8) < Tc_min
                mod_com(8) = Tc_min;
            end

            if mod_com(1)>A_max 
                mod_com(1) = A_max;
            end
            if mod_com(2)> B_max
                mod_com(2) = B_max;
            end
            if mod_com(3)> C_max
                mod_com(3) = C_max;
            end
            if mod_com(4)>D_max
                mod_com(4) = D_max;
            end
            if mod_com(5)>V0_max
                mod_com(5) = V0_max;
            end
            if mod_com(6)>Ta_max
                mod_com(6) = Ta_max;
            end
            if mod_com(7) > Tb_max
                mod_com(7) = Tb_max;
            end
            if mod_com(8) > Tc_max
                mod_com(8) = Tc_max;
            end
        %Perhitungan misfit untuk prosedur komensalisme
        [cal_slip_com] = cal_slip(mod_com(1),mod_com(2),mod_com(3),mod_com(4),mod_com(5),mod_com(6),mod_com(7),mod_com(8),t);
         err_com = fun_obj(Us_data,cal_slip_com);
         %Update model jika  nilai misfit lebih baik proses komensalisme
         if err_com < E(i)
             model(i,:) = mod_com;
             E(i) = err_com;
             Us_cal(i,:) = cal_slip_com;
         end
         %Parasitisme
         j = randi(npop,1);
         if j == i 
             j = randi(npop,1);
        end
         mod_par = model(i,:);
         p1 = randi(8,1);
         if p1 == 1
            mod_par(i,1) = A_min + rand*(A_max-A_min);
         elseif p1 == 2
            mod_par(i,2) = B_min + rand*(B_max-B_min);
         elseif p1 == 3
             mod_par(i,3) = C_min + rand*(C_max-C_min);
         elseif p1 == 4
             mod_par(i,4) = D_min + rand*(D_max-D_min);
         elseif p1 == 5
             mod_par(i,5) = V0_min + rand*(V0_max-V0_min);
         elseif p1 == 6
             mod_par(i,6) = Ta_min + rand*(Ta_max-Ta_min);
         elseif p1 == 7
             mod_par(i,7) = Tb_min + rand*(Tb_max-Tb_min);
         elseif p1 == 8
             mod_par(i,8) = Tc_min + rand*(Tc_max-Tc_min);
         end
         %Perhitungan misfit untuk tahap parasitisme
         [cal_slip_par] = cal_slip(mod_par(1),mod_par(2),mod_par(3),mod_par(4),mod_par(5),mod_par(6),mod_par(7),mod_par(8),t);
         err_par = fun_obj(Us_data,cal_slip_par);
         %Update model jika  nilai misfit lebih baik proses komensalisme
         if err_par < E(i)
             model(j,:) = mod_par(1,:);
             E(j) = err_par;
             Us_cal(i,:) = cal_slip_par;
         end
    end
    %Update model terbaik untuk setiap iterasi
    Emin = 100;
    for ipop = 1 : npop
        Emin = E(ipop);
        model_baru = model(ipop,:);
        Us_model = Us_cal(ipop,:);
    end
    %Nilai misfit terbaik
    Egen(itr)=Emin;
    % Plotting 
figure(1)
scatter(t,Us_data,'b.')
hold on
plot(t,Us_model,'r')
xlabel('Day','FontSize',10,'FontWeight','Bold')
ylabel('Displacement Us (cm)','FontSize',10,'FontWeight','Bold')
title(['Exponential and Logarithm models of Easting LEWK|| Hasil Inversi ke-.',num2str(itr)],'FontWeight','bold')
subtitle(['A = ',num2str(model_baru(1)),' ; B = ',num2str(model_baru(2)),' ; C = ',num2str(model_baru(3)),' ; D = ',num2str(model_baru(4)),' ; V0 = ',num2str(model_baru(5)),' ; Ta = ',num2str(model_baru(6)),' ; Tb = ',num2str(model_baru(7)),' ; Tc = ',num2str(model_baru(8)), '|| ERMS = ',num2str(Egen(itr))],'FontWeight','bold')
legend('GPS data','models')
grid on
set(gcf, 'Position', get(0, 'Screensize'));
hold off
frame = getframe(gcf);
writeVideo(v,frame);
end
close(v);
%Plot grafik misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Misfit Graph ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on

saveas(figure(2),'Grafik misfit SOS.png')
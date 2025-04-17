%% Получение и визуализация данных лабораторного эксперимента
eRange = -0.49:0.01:0.98;
sigm = [0.633 0.4521 0.2713 0.0904];

try
    x_step = 17.65;
    y_step = 405.6;
    x_shift = [40.070745 40.070745 41.253374 41.844688];
    y_shift = [4.1597764 4.1597764 5.7850272 4.1597764];
    r = 11.6; 
    eps_init = (y_shift./y_step)./r;

    black = (importdata('C:\YandexDisk\Science\icmmLaboratory\Graph2Digit\graph_1_data\Points\black.dat'))';
    blue = (importdata('C:\YandexDisk\Science\icmmLaboratory\Graph2Digit\graph_1_data\Points\blue.dat'))';
    green = (importdata('C:\YandexDisk\Science\icmmLaboratory\Graph2Digit\graph_1_data\Points\green.dat'))';
    red = (importdata('C:\YandexDisk\Science\icmmLaboratory\Graph2Digit\graph_1_data\Points\red2.dat'))';
    %black: 0.633 Pa, blue: 0.4521 Pa, green: 0.2713 Pa, red: 0.0904 Pa
    x_black = ((black(1,:)-x_shift(1)) / x_step)';
    y_black = (((black(2,:)) / y_step) / r)';
    x_blue = ((blue(1,:)-x_shift(2)) / x_step)';
    y_blue = (((blue(2,:)) / y_step) / r)';
    x_green = ((green(1,:)-x_shift(3)) / x_step)';
    y_green = (((green(2,:)) / y_step) / r)';
    x_red = ((red(1,:)-x_shift(4)) / x_step)';
    y_red = (((red(2,:)) / y_step) / r)';

    p = plot(x_black,y_black,'k',x_blue,y_blue,'b', x_green,y_green,'g',x_red,y_red,'r');
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p(1).LineWidth = 2;
    p(2).LineWidth = 2;
    p(3).LineWidth = 2;
    p(4).LineWidth = 2;
    ylim([0 0.11])
    legend({'\sigma = 0.633 Pa','\sigma = 0.4521 Pa', '\sigma = 0.2713 Pa', '\sigma = 0.0904 Pa'},'Location','northwest');
 
    x_black = x_black([1 6:end]);
    y_black = y_black([1 6:end]);
    x_blue = x_blue([1 4:end]);
    y_blue = y_blue([1 4:end]);
    x_green = x_green([1 4:end]);
    y_green = y_green([1 4:end]);
    x_red = x_red([1 4:end]);
    y_red = y_red([1 4:end]);
    time = {x_black, x_blue, x_green, x_red};
    excpt = false;
catch
    disp("Cannot find experimental data. Make sure that file names is correct");
    eps_init = [0.001 0.001 0.001 0.001];
    time = {0:0.05:8; 0:0.05:8; 0:0.05:8; 0:0.05:8};
    excpt = true;
end

%% Идентификация параметров математической модели для напряжения sigm = 0.633 Па

t = time{1};
thet = 300 * 1.38 * 10^-23;
lamd = 200e-19;
G0 = 170;
G1 = 20.6e6;
tau0 = 30.13e-10;
gam = 0.5e-21; %начальные значения параметров

eps = FindDeformation(t, G0, G1, lamd, gam, tau0, thet, sigm(1), eps_init(1), eRange); %определение зависимости eps(t) с нач. значениями

if excpt %если нет экспериментальных данных
    initPar = [G0 G1 lamd gam tau0];
    figure();
    p1 = plot(t, eps,'k'); %сопоставление численных (при начальных параметрах) и экспериментальных результатов
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p1(1).LineWidth = 1;
    ylim([0 inf]);
    xlim([0 inf]);
    grid on;
    legend({'Numerical results with initial parameters for \sigma = 0.633 Pa'},'Location','southeast');
else
    figure();
    p1 = plot(x_black,y_black,'ko',t, eps,'k'); %сопоставление численных (при начальных параметрах) и экспериментальных результатов
        xlabel('t', FontSize=14);
        ylabel('\epsilon(t)', FontSize=14);
        p1(1).LineWidth = 0.7;
        p1(1).MarkerSize = 1.8;
        p1(2).LineWidth = 2;
        ylim([0 inf]);
        xlim([0 inf]);
        grid on;
        legend({'Experimental data for \sigma = 0.633 Pa','Numerical results with initial parameters'},'Location','southeast');
    
    x0 = [G0 G1 lamd gam tau0]; %решение оптимизационной задачи и идентификация параметров модели
    options = optimset('MaxIter', 10^4,'MaxFunEvals',10^4);
    [xSol_black, valOF] = fminsearch(@(x) FindObjectiveFunction(x_black, y_black, x(1), x(2), x(3), x(4), x(5), thet, sigm(1), eps_init(1),eRange), x0, options);
    eNum_black = FindDeformation(t, xSol_black(1), xSol_black(2), xSol_black(3), xSol_black(4), xSol_black(5), thet, sigm(1), eps_init(1),eRange);
    initPar = xSol_black;

    figure();
    p2 = plot(x_black,y_black,'ko',t, eNum_black,'k'); %сопоставление численных (при идентифицированных параметрах) и экспериментальных результатов
        xlabel('t', FontSize=14);
        ylabel('\epsilon(t)', FontSize=14);
        p2(1).LineWidth = 0.7;
        p2(1).MarkerSize = 1.8;
        p2(2).LineWidth = 2;
        grid on;
        legend({'Experimental data for \sigma = 0.633 Pa','Numerical results with identified parameters'},'Location','southeast');
end

%% Идентификация параметров математической модели для напряжения sigm = 0.4521 Па

t = time{2};
lamd = initPar(3);
G0 = initPar(1);
G1 = initPar(2);
tau0 =  initPar(5);
gam = initPar(4);

eps = FindDeformation(t, G0, G1, lamd, gam, tau0, thet, sigm(2), eps_init(2),eRange); 

if excpt %если нет экспериментальных данных
    figure();
    p1 = plot(t, eps,'b'); %сопоставление численных (при начальных параметрах) и экспериментальных результатов
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p1(1).LineWidth = 1;
    ylim([0 inf]);
    xlim([0 inf]);
    grid on;
    legend({'Numerical results with initial parameters for \sigma = 0.4521 Pa'},'Location','southeast');
else
    figure();
    p1 = plot(x_blue,y_blue,'bo',t, eps,'b'); 
        xlabel('t', FontSize=14);
        ylabel('\epsilon(t)', FontSize=14);
        p1(1).LineWidth = 0.7;
        p1(1).MarkerSize = 1.8;
        p1(2).LineWidth = 2;
        ylim([-inf inf])
        xlim ([-inf inf])
        grid on;
        legend({'Experimental data for \sigma = 0.4521 Pa','Numerical results with initial parameters'},'Location','southeast');
    
    x0 = [G0 G1 lamd gam tau0]; 
    options = optimset('MaxIter', 10^4,'MaxFunEvals',10^4);
    [xSol_blue, valOF] = fminsearch(@(x) FindObjectiveFunction(x_blue, y_blue, x(1), x(2), x(3), x(4), x(5), thet, sigm(2), eps_init(2),eRange), x0, options);
    eNum_blue = FindDeformation(t, xSol_blue(1), xSol_blue(2), xSol_blue(3), xSol_blue(4), xSol_blue(5), thet, sigm(2), eps_init(2),eRange); 
    
    figure();
    p2 = plot(x_blue,y_blue,'bo',t, eNum_blue,'b'); 
        xlabel('t', FontSize=14);
        ylabel('\epsilon(t)', FontSize=14);
        p2(1).LineWidth = 0.7;
        p2(1).MarkerSize = 1.8;
        p2(2).LineWidth = 2;
        grid on;
        legend({'Experimental data for \sigma = 0.4521 Pa','Numerical results with identified parameters'},'Location','southeast');
end

%% Идентификация параметров математической модели для напряжения sigm = 0.2713 Па

t = time{3};
lamd = initPar(3);
G0 = initPar(1);
G1 = initPar(2);
tau0 =  initPar(5);
gam = initPar(4);

eps = FindDeformation(t, G0, G1, lamd, gam, tau0, thet, sigm(3), eps_init(3),eRange); 

if excpt %если нет экспериментальных данных
    figure();
    p1 = plot(t, eps,'g'); %сопоставление численных (при начальных параметрах) и экспериментальных результатов
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p1(1).LineWidth = 1;
    ylim([0 inf]);
    xlim([0 inf]);
    grid on;
    legend({'Numerical results with initial parameters for \sigma = 0.2713 Pa'},'Location','southeast');
else
    figure();
    p1 = plot(x_green,y_green,'go',t, eps,'g'); 
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p1(1).LineWidth = 0.7;
    p1(1).MarkerSize = 1.8;
    p1(2).LineWidth = 2;
    ylim([0 inf])
    grid on;
    legend({'Experimental data for \sigma = 0.2713 Pa','Numerical results with initial parameters'},'Location','southeast');

    x0 = [G0 G1 lamd gam tau0];
    options = optimset('MaxIter', 10^4,'MaxFunEvals',10^4);
    [xSol_green, valOF] = fminsearch(@(x) FindObjectiveFunction(x_green, y_green, x(1), x(2), x(3), x(4), x(5), thet, sigm(3), eps_init(3),eRange), x0, options);
    eNum_green = FindDeformation(t, xSol_green(1), xSol_green(2), xSol_green(3), xSol_green(4), xSol_green(5), thet, sigm(3), eps_init(3),eRange); 
    
    figure();
    p2 = plot(x_green,y_green,'go',t, eNum_green,'g'); 
        xlabel('t', FontSize=14);
        ylabel('\epsilon(t)', FontSize=14);
        p2(1).LineWidth = 0.7;
        p2(1).MarkerSize = 1.8;
        p2(2).LineWidth = 2;
        grid on;
        legend({'Experimental data for \sigma = 0.2713 Pa','Numerical results with identified parameters'},'Location','southeast');
end

%% Идентификация параметров математической модели для напряжения sigm = 0.0904 Pa

t = time{4};
lamd = initPar(3);
G0 = initPar(1);
G1 = initPar(2);
tau0 =  initPar(5);
gam = initPar(4);

eps = FindDeformation(t, G0, G1, lamd, gam, tau0, thet, sigm(4), eps_init(4),eRange); 

if excpt %если нет экспериментальных данных
    figure();
    p1 = plot(t, eps,'r'); %сопоставление численных (при начальных параметрах) и экспериментальных результатов
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p1(1).LineWidth = 1;
    ylim([0 inf]);
    xlim([0 inf]);
    grid on;
    legend({'Numerical results with initial parameters for \sigma = 0.0904 Pa'},'Location','southeast');
else
    figure();
    p1 = plot(x_red,y_red,'ro',t, eps,'r'); 
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p1(1).LineWidth = 0.7;
    p1(1).MarkerSize = 1.8;
    p1(2).LineWidth = 2;
    grid on;
    ylim([0 inf])
    legend({'Experimental data for \sigma = 0.0904 Pa','Numerical results with initial parameters'},'Location','southeast');
    
    x0 = [G0 G1 lamd gam tau0]; 
    options = optimset('MaxIter', 10^4,'MaxFunEvals',10^4);
    [xSol_red, valOF] = fminsearch(@(x) FindObjectiveFunction(x_red, y_red, x(1), x(2), x(3), x(4), x(5), thet, sigm(4), eps_init(4),eRange), x0, options);
    eNum_red = FindDeformation(t, xSol_red(1), xSol_red(2), xSol_red(3), xSol_red(4), xSol_red(5), thet, sigm(4), eps_init(4),eRange); 

    figure();
    p2 = plot(x_red,y_red,'ro',t, eNum_red,'r'); 
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p2(1).LineWidth = 0.7;
    p2(1).MarkerSize = 1.8;
    p2(2).LineWidth = 2;
    grid on;
    legend({'Experimental data for \sigma = 0.0904 Pa','Numerical results with identified parameters'},'Location','southeast');    
end 

%% Вычисление среднего и среднеквадратичного отклонения параметров модели, построение численных и экспериментальных результатов

if (excpt == false) 
    identPar = [xSol_black; xSol_blue; xSol_green; xSol_red]';
    meanPar = [mean(identPar(1,:)) mean(identPar(2,:)) mean(identPar(3,:)) mean(identPar(4,:)) mean(identPar(5,:))];
    stdDevPar = [std(identPar(1,:)) std(identPar(2,:)) std(identPar(3,:)) std(identPar(4,:)) std(identPar(5,:))];
    
    p = plot(x_black,y_black,'ko',x_black, eNum_black,'r', x_blue,y_blue,'ko',  x_blue, eNum_blue,'r', ...
        x_green,y_green,'ko',x_green, eNum_green,'r', x_red, y_red,'ko',x_red, eNum_red,'r');
        xlabel('t, c', FontSize=12);
        ylabel(char(949));
        grid on;
        ylim([0 0.105]);
        xlim([0 8]);
    p(1).MarkerSize = 3;
    p(1).MarkerFaceColor = 'k';
    p(3).MarkerSize = 3;
    p(3).MarkerFaceColor = 'k';
    p(5).MarkerSize = 3;
    p(5).MarkerFaceColor = 'k';
    p(7).MarkerSize = 3;
    p(7).MarkerFaceColor = 'k';
    p(2).LineWidth = 2.0;
    p(4).LineWidth = 2.0;
    p(6).LineWidth = 2.0;
    p(8).LineWidth = 2.0;
    
    %вычисление относительной погрешности для каждой кривой
    w_black = max((abs(y_black - eNum_black)./y_black) * 100);
    w_blue = max((abs(y_blue - eNum_blue)./y_blue) * 100);
    w_green = max((abs(y_green - eNum_green)./y_green) * 100);
    w_red = max((abs(y_red - eNum_red)./y_red) * 100);
end

%% Построение дополнительных графиков при sigm = 0.633 Па

if (excpt == false)
    sigmVal = 0.633; 
    thet = 300 * 1.38 * 10^-23;
    G0 = xSol_black(1);
    G1 = xSol_black(2);
    lamd = xSol_black(3);
    gam = xSol_black(4);
    tau0 = xSol_black(5);
    chi = thet / (lamd * gam);
    t = 0:0.05:8;
    e_init = eps_init(1);
    eTotal = FindDeformation(t, G0, G1, lamd, gam, tau0, thet, sigmVal, e_init,eRange); 
    eOrient = eTotal - ((sigmVal - G0 * eTotal) / (G1)); %деформация epso, обусловленная ориентационными эффектами
    sOrient = sigmVal - eTotal * G0; %напряжения ориентационного элемента в схеме модели
    
    % построение зависимости составляющих полной деформации от времени
    figure();
    p1 = plot(t, eTotal, 'b', t, eOrient, '--r', t, eTotal-eOrient); 
        xlabel('t', FontSize=14);
        ylabel('Deformation', FontSize=14);
        p1(1).LineWidth = 2;
        p1(2).LineWidth = 3.5;
        p1(3).LineWidth = 2;
        grid on;
        ylim([0 0.105]);
        xlim([0 8]);
        legend({char(949), [char(949) '_{o}'], [char(949) '_{e}']}, 'FontSize', 11, 'Location','southeast');
    
    % построение зависимости составляющих полных напряжений от времени
    sTotal = sigmVal * ones(1,length(t));
    figure();
    p2 = plot(t, sTotal, 'b', t, sOrient, 'r', t, sTotal - sOrient);
        xlabel('t', FontSize=14);
        ylabel('Stress, Pa', FontSize=14);
        p2(1).LineWidth = 2;
        p2(2).LineWidth = 2;
        p2(3).LineWidth = 2;
        grid on;
        ylim([-inf inf]);
        xlim([0 inf]);
        legend({char(963), [char(963) '_{o}'], [char(963) '_{r}']}, 'Location','northeast');
    
    % построение зависимости потенциала Psi от параметров epso и sigmo
    tCase = [0; 0.4; 5];
    FCase = [Fpotential(sOrient(t==tCase(1)) * (gam/thet), eRange, chi);
             Fpotential(sOrient(t==tCase(2)) * (gam/thet), eRange, chi);
             Fpotential(sOrient(t==tCase(3)) * (gam/thet), eRange, chi)];
    eOCase = [eOrient(t==tCase(1)) eOrient(t==tCase(2)) eOrient(t==tCase(3))];
    FCasePoint = [[0 0.05 0.095]; [FCase(1,50) FCase(2, 55) FCase(3,59)]];
    figure();
    p3 = plot(eRange,FCase(1,:), eRange, FCase(2,:),  eRange,FCase(3,:), ...
        FCasePoint(1,1), FCasePoint(2,1), 'o', FCasePoint(1,2), FCasePoint(2,2), 'o', FCasePoint(1,3), FCasePoint(2,3), 'o');
        xlim([min(eRange) max(eRange)])
        ylim([-inf inf])
        xlabel([char(949) '_{o}'], FontSize=14);
        ylabel('\Psi-<\Psi>', FontSize=14);
        p3(1).LineWidth = 2;
        p3(2).LineWidth = 2;
        p3(3).LineWidth = 2;
        p3(4).LineWidth = 2;
        p3(5).LineWidth = 2;
        p3(6).LineWidth = 2;
        grid on;
        ylim([-inf inf]);
        xlim([-inf inf]);
        legend({['\Psi(' char(963) '_{o}(0),' char(949) '_{o}(t)'], ['\Psi(' char(963) '_{o}(0.4),' char(949) '_{o}(t)'], ...
           ['\Psi(' char(963) '_{o}(5),' char(949) '_{o}(t)'], ['\Psi(' char(963) '_{o}(0),' char(949) '_{o}(0)'], ...
          ['\Psi(' char(963) '_{o}(0.4),' char(949) '_{o}(0.4)'], ['\Psi(' char(963) '_{o}(5),' char(949) '_{o}(5)']}, 'Location','northeast');
    
    % построение зависимости величины барьера от времени нагружения при
    % различных значениях внешних напряжений
    nu = [];
    dF = [];
    dFCur = [];
    nuCur = []; 
    thet = 300 * 1.38 * 10^-23;
    G0 = meanPar(1);
    G1 = meanPar(2);
    lamd = meanPar(3);
    gam = meanPar(4);
    tau0 = meanPar(5);
    sigmVal = [0.633 0.4521 0.2713 0.0904];
    e_init = mean(eps_init);
    
    for k=1:length(sigmVal)
        eTotal = FindDeformation(t, G0, G1, lamd, gam, tau0, thet, sigmVal(k), e_init,eRange);
        sOrient = sigmVal(k) - eTotal * G0;
        for i=1:length(sOrient)
            sOrientCur = sOrient(i);
            dFCur(end + 1) = deltaF(sOrientCur * (gam/thet), eRange, chi);
            nuCur(end + 1) = nuCalculation(sOrientCur * (gam/thet), eRange, chi, G1, tau0);
        end
        dF{end + 1} = dFCur;
        nu{end + 1} = nuCur;
        dFCur = [];
        nuCur = [];
    end
    
    f1 = figure();
    p4 = plot(t, dF{1}, '--', t, dF{2}, t, dF{3}, t, dF{4},'--'); 
        xlabel('t, s', FontSize=14);
        ylabel('\Delta\Psi', FontSize=14);
        p4(1).LineWidth = 1.5;
        p4(2).LineWidth = 1.5;
        p4(3).LineWidth = 1.5;
        p4(4).LineWidth = 1.5;
        grid on;
        legend({[char(963) ' = 0.633 Pa'],[char(963) ' = 0.4521 Pa'],[char(963) ' = 0.2713 Pa'],[char(963) ' = 0.0904 Pa']});
    axes('position', [0.35 0.25 0.3 0.3]);
    box on
    your_index = 0.2<t & t<0.5;
    plot(t(your_index), dF{1}(your_index), '--', t(your_index), dF{2}(your_index), ...
        t(your_index), dF{3}(your_index),t(your_index), dF{4}(your_index),'--');
    axis tight
end

%% локальные функции
function valOF = FindObjectiveFunction(tExp, eExp, G0, G1, lamd, gam, tau0, thet, sigm, epsInit,eRange) %найти значение целевой функции при известных экспериментальных деформациях eExp и векторе параметров модели xParam
eNum = FindDeformation(tExp, G0, G1, lamd, gam, tau0, thet, sigm, epsInit,eRange);
valOF = sqrt(sum(((eExp-eNum)./eExp).^2));
end

function def = FindDeformation(t, G0, G1, lamd, gam, tau0, thet, sigm, epsInit, epsRange) %зависимость общей деформации от времени при постоянном напряжении
[~, eps] = ode15s(@(t, eps) odeF(eps, G0, G1, lamd, gam, tau0, thet, sigm, epsRange), t, epsInit);
def = eps;
end

function dedt = odeF(eps, G0, G1, lamd, gam, tau0, thet, sigm, epsRange) %функция производной для решения ДУ в FindDeformation
chi = thet / (lamd * gam);
eps0 = eps - ((sigm - G0 * eps) / (G1));
sigm0 = sigm - eps * G0;
ksiVal = ksiDependence(eps0);
nu = nuCalculation(sigm0 * (gam/thet),epsRange, chi, G1, tau0);
dedt = (G1 / (nu * (G1 + G0))) * (sigm - G0 * eps - (thet/(gam)) * ksiVal + (lamd/gam) * eps0);
end

function ksi = ksiDependence(eps0)
    ksi = (260 + eps0 .* (123200 + (76250-114200 * eps0) .* eps0))./(12400 + eps0 .* (16770 + eps0 .* (-20670 + (-8820 + eps0) .* eps0)));
end

function Fp = Fpotential(s, eps, chi) %вычисление потенциала Psi 
    Fp = real(-1.7*log(0.99-eps) - 114161*log(8822.34-eps) - 0.9*log(0.5+eps) - 36.34*log(2.84+eps) - (eps.^2)/(2*chi) - s.* eps);
    means = mean(Fp, 2);
    Fp = Fp - means;
end

function dF = deltaF(s, eps, chi) %вычисление величины барьера потенциала Psi
F = Fpotential(s,eps,chi);

A1 = islocalmax(F);
A2 = islocalmin(F);

if(numel(eps(A1)) + numel(eps(A2)) == 0)
    if(F(1) > F(end))
        dF = 0;
    else
        dF = F(end) - F(1);
    end
else
    if(numel(eps(A1)) + numel(eps(A2)) > 1) 
        MaxInd = find(A1, 1); 
        MinInd = find(A2, 1); 
        if(eps(MaxInd) > eps(MinInd))
            dF = F(MaxInd) - F(MinInd);
        else 
            disp("Error: Xmax < Xmin");
        end
    else
        if(F(1) > F(end))
            dF = 0;
        else
            dF = F(end) - F(A2);
        end
    end
end
end

function nu = nuCalculation(s0d, e, chi, G1, tau0) % вычисление ориентационной вязкости
    nu = G1 * tau0 * exp(deltaF(s0d, e, chi));
end

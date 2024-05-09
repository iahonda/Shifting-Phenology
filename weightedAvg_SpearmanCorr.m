%% Statistical Analysis for manuscript "Shifting Phenology as a Key Driver of Shelf Zooplankton Population Variability" (submitted to Ecology)
% MATLAB code for calculating weighted averages and plotting figures
% Isabel Honda
% May 2024

%% Calculating Weighted Average Curves and Spearman Rank Correlation Plots

clear; clc; close all

var1_all = % Load in data for variable 1 (e.g., C. finmarchicus abundance)
var2_all = % Load in data for variable 2 (e.g., C. finmarchicus abundance)
% Putting both variables on the same year scale
yearValues = intersect(unique(var1_all.year),unique(var2_all.year));
var1 = var1_all(var1_all.year >= yearValues(1) & var1_all.year <= yearValues(end),:);
var2 = var2_all(var2_all.year >= yearValues(1) & var2_all.year <= yearValues(end),:);

EcoMon = % Load in MARMAP/EcoMon raw data
EcoMon_basins = EcoMon(EcoMon.strata_47 == 37 | EcoMon.strata_47 == 42 | EcoMon.strata_47 == 38 | EcoMon.strata_47 == 39,["strata_47","year","doy","calfin_100m3"]);
logEcoMon_basins = EcoMon_basins;
logEcoMon_basins.calfin_100m3 = log(EcoMon_basins.calfin_100m3/100 + 1);
logEcoMon_basins = rmmissing(logEcoMon_basins);

% For plotting
var1_name = % Varible 1 name (e.g., C. finmarchicus Spring abundance)
var2_name = % Varible 2 name (e.g., C. finmarchicus Fall abundance)
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]}; 
c1 = colors{1};
c2 = colors{2};
colorsR = {'#6d636b','#77926f','#b67057','#c86328','#b2af28'};
cR=1;

% Winter: 1 to 90 | Spring: 91 to 181 | Summer: 182 to 272 | Fall: 273 to 365
filter_var1 = % Choose variable 1 filter (e.g., Spring DoYs - var1.doy >= 91 & var1.doy <= 181;)
filter_var2 = % Choose variable 2 filter (e.g., Fall DoYs - var2.doy >= 274 & var2.doy <= 365;) 
strataNum = {37, 42, 48}; % For WB, JB, and GBn, respectively
strata = {'s37','s42','s48'};
ses = {'se_37','se_42','se_48'};
var1s = table();
var2s = table();

% Weighted average function
weighted_avg = @(x, se) sum(x .* (1 ./ se.^2)) / sum(1 ./ se.^2);

titles = {'Wilkinson Basin','Jordan Basin','Georges Basin'};
figure(1)
clf
set(gcf,'color','w');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

clear h
springs = table();
falls = table();

for i=1:length(strata)
    
    mean_var2 = splitapply(weighted_avg, var2.(strata{i})(filter_var2), var2.(ses{i})(filter_var2),findgroups(var2.year(filter_var2)));
    mean_var1 = splitapply(weighted_avg, var1.(strata{i})(filter_var1), var1.(ses{i})(filter_var1),findgroups(var1.year(filter_var1)));

    % Saving means for outside of loop
    var2s.(strata{i}) = mean_var2;
    var1s.(strata{i}) = mean_var1;

    
    subplot(3,2,2*i-1);
    yyaxis left
    hold on
    h(1) = plot(yearValues,mean_var1,'LineWidth',2,'Color',c1)
    ylabel(var1_name)
    ax = gca;
    ax.YAxis(1).Color = c1;

    yyaxis right
    h(2) = plot(yearValues,mean_var2,'LineWidth',2,'Color',c2)
    ylabel(var2_name)
    ax = gca;
    ax.YAxis(2).Color = c2;

    minSpring = 10;
    maxSpring = 0;
    minFall = 10;
    maxFall = 0;
    
    for y=yearValues(1):yearValues(end)
        if strataNum{i} == 48
            spring = logEcoMon_basins(logEcoMon_basins.year == y & logEcoMon_basins.doy >= 91 & logEcoMon_basins.doy <= 181 & (logEcoMon_basins.strata_47 == 38 | logEcoMon_basins.strata_47 == 39),"calfin_100m3");    
            fall = logEcoMon_basins(logEcoMon_basins.year == y & logEcoMon_basins.doy >= 274 & logEcoMon_basins.doy <= 365 & (logEcoMon_basins.strata_47 == 38 | logEcoMon_basins.strata_47 == 39),"calfin_100m3");
            spring_doys = logEcoMon_basins(logEcoMon_basins.year == y & logEcoMon_basins.doy >= 91 & logEcoMon_basins.doy <= 181 & (logEcoMon_basins.strata_47 == 38 | logEcoMon_basins.strata_47 == 39),"doy");    
            fall_doys = logEcoMon_basins(logEcoMon_basins.year == y & logEcoMon_basins.doy >= 274 & logEcoMon_basins.doy <= 365 & (logEcoMon_basins.strata_47 == 38 | logEcoMon_basins.strata_47 == 39),"doy");
        else
            spring = logEcoMon_basins(logEcoMon_basins.year == y & logEcoMon_basins.doy >= 91 & logEcoMon_basins.doy <= 181 & logEcoMon_basins.strata_47 == strataNum{i},"calfin_100m3");    
            fall = logEcoMon_basins(logEcoMon_basins.year == y & logEcoMon_basins.doy >= 274 & logEcoMon_basins.doy <= 365 & logEcoMon_basins.strata_47 == strataNum{i},"calfin_100m3");
            spring_doys = logEcoMon_basins(logEcoMon_basins.year == y & logEcoMon_basins.doy >= 91 & logEcoMon_basins.doy <= 181 & logEcoMon_basins.strata_47 == strataNum{i},"doy"); 
            fall_doys = logEcoMon_basins(logEcoMon_basins.year == y & logEcoMon_basins.doy >= 274 & logEcoMon_basins.doy <= 365 & logEcoMon_basins.strata_47 == strataNum{i},"doy");
        end

        minSpringc = min(spring.calfin_100m3(spring.calfin_100m3>0));
        if minSpringc < minSpring
            minSpring = minSpringc;
        end
        maxSpringc = max(spring.calfin_100m3);
        if maxSpringc > maxSpring
            maxSpring = maxSpringc;
        end
    
        minFallc = min(fall.calfin_100m3(fall.calfin_100m3>0));
        if minFallc < minFall
            minFall = minFallc;
        end
        maxFallc = max(fall.calfin_100m3);
        if maxFallc > maxFall
            maxFall = maxFallc;
        end

    if isempty(spring)
        continue
    else
        springs.(strata{i})(y-1976) = mean(spring.calfin_100m3);
        yyaxis left
        plot(y+spring_doys.doy/365,(spring.calfin_100m3),'.','color',[0 0.4470 0.7410]);
        plot(y+5/12,mean(spring.calfin_100m3),'.','markersize',20,'color',[0 0.4470 0.7410]);
    end
    if isempty(fall)
        continue
    else
        yyaxis right
        falls.(strata{i})(y-1976) = mean(fall.calfin_100m3);
        plot(y+fall_doys.doy/365,(fall.calfin_100m3),'.','color',[0.8500 0.3250 0.0980]);
        plot(y+11/12,mean(fall.calfin_100m3),'.','markersize',20,'color',[0.8500 0.3250 0.0980]);
        hold on
    end

    end
    
    h(1) = plot(NaN,NaN,'k.');
    h(2) = plot(NaN,NaN,'k.','MarkerSize',20);
    h(3) = plot(NaN,NaN,'k-','LineWidth',2);

    yyaxis left
    ylim([minSpring maxSpring])

    yyaxis right
    ylim([minFall maxFall])

    xlim([min(yearValues) max(yearValues)])
    xlabel('Year')
    title(titles{i})
    if i==1
        legend(h,'MARMAP/EcoMon Data','MARMAP/EcoMon Average','GAM Average','location','se')
    end
    ax = gca; 
    ax.FontSize = 15;
end


for i=1:length(strata)
    X = var1s.(strata{i});
    Y = var2s.(strata{i});

    sF = [X Y];
    sF = sF(all(sF,2),:);
    
    % Compute orthogonal regression slope (weighted)
    [rho,pval] = corr(sF(:,1),sF(:,2),'type','Spearman');
  
    % Plotting the orthogonal regression line with error bars
    a(i) = subplot(3,2,2*i)
    scatter(sF(:,1),sF(:,2),[],'MarkerFaceColor',colorsR{cR},'MarkerEdgeColor',colorsR{cR},'Marker','square')
    if pval < 0.001
        txt = {['\rho = ' num2str(rho, '%.4f')],'{\it p} value < 0.001'};
    else
        txt = {['\rho = ' num2str(rho, '%.4f')],['{\it p} value = ' num2str(pval, '%.4f')]};
    end
    annotation('textbox','String',txt,'Position',a(i).Position,'Vert','bottom','FitBoxToText','on', 'Color', 'r','FontSize', 15)
    xlabel(var1_name)
    ylabel(var2_name)
    title(titles{i})
    grid on;
    hold off;
    ax = gca; 
    ax.FontSize = 15;
end
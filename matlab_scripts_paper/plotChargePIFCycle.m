clear all;

colormap default;
set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex'); 

nranks = 16;
sranks = 4;
nCycles = 1;
time_str = 'T_192_dt_003125'
grid_str = '64_cube';
Np_str =  'Pc_10';
test_str = 'TSI';
dir = ['../../', test_str,'/corrected_shape_function/Conservation_studies/',time_str,'/', Np_str, '/', grid_str, '/'];
color_map = get(0, 'DefaultAxesColorOrder');

if(strcmp(test_str,'LandauDamping') || strcmp(test_str,'TSI'))
    initial_charge = -(4*pi)^3;
else
    initial_charge = -1562.5;
end
initial_charge
iterRank = cell(nranks, nCycles);

fig=figure;
for nc=1:nCycles
    for r=1:nranks
        file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/para_tol_1em8/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
        B = readmatrix(file,'NumHeaderLines',1,'Delimiter',' ');
        iterRank{r,nc} = B(:,1);
    end
end

[val, idx] = max(max(cellfun(@max, iterRank)));


for nc=1:nCycles
    if(mod(nc,2) == 0)
        lastRank = 1;
        start_rank = nranks;
        end_rank = 1;
        step = -1;
    else
        lastRank = nranks;
        start_rank = 1;
        end_rank = nranks;
        step = 1;
    end
    max_iter = iterRank{lastRank,nc}(end);
    
    chargeIter = cell(max_iter, 1);
    timeIter   = cell(max_iter, 1);
    for iter=1:max_iter
        shift=1;
        total = 0;
        timeIter{iter} = zeros(1,1);
        chargeIter{iter} = zeros(1,1);
        for r=start_rank:step:end_rank
            r
            if(iter <= iterRank{r,nc}(end))
                file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/para_tol_1em8/data/Energy_rank_', num2str(r-1),'_nc_',num2str(nc),'_iter_',num2str(iter),'.csv'];
                B = readmatrix(file,'NumHeaderLines',0,'Delimiter',' ');
                total = total + size(B,1);
                chargeIter{iter}(shift:total, :) = B(:,5);
                timeIter{iter}(shift:total, :) = B(:,1);
                shift = shift + size(B,1);
            end
        end
    end
    
    
    for iter=1:max_iter
        charge_error_pif = abs(chargeIter{iter}(:,1) - initial_charge)./abs(initial_charge);
    
        if(iter <= 7)
            nco = iter;
            linestyle = '-';
        elseif((iter > 7) && (iter <=14))
            nco = mod(iter, 7);
            if(nco == 0)
                nco = 1;
            end
            linestyle = '--';
        elseif((iter > 14) && (iter <=21))
            nco = mod(iter, 7);
            if(nco == 0)
                nco = 1;
            end
            linestyle = '-.';
        end
    
        pl(nc, iter) = semilogy(timeIter{iter}(:),charge_error_pif(:),'LineStyle',linestyle,'Color',color_map(nco,:),'LineWidth',1.5);
        hold on;
    end
end

%%Plot reference charge error curve for serial time stepping
dir_serial = ['../../../ElectrostaticPIF/',test_str,'_conservation_studies/corrected_shape_function/'];
A_pif=readmatrix([dir_serial,'64_64_64_Pc_10/T_192/ngpus_',num2str(sranks),...
                  '/dt_003125/fine_tol_1em7/data/Energy_',num2str(sranks),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
charge_error_pif_serial = abs(A_pif(:,5) - initial_charge)/abs(A_pif(1,5));
semilogy(A_pif(:,1),charge_error_pif_serial(:),'k--','LineWidth',2.0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hold off;
grid on;
xlabel('time');
%ylabel('Rel. Charge error');
set(gca,'Fontsize',22);
%legend(pl(idx,:),'k = 1','k = 2','k = 3','k = 4','k = 5','k = 6','k = 7','k = 8',...
%                 'k = 9','k = 10','k = 11','k = 12','k = 13','k = 14','k = 15','k = 16',...
%                 'Location','northoutside','Orientation','horizontal','NumColumns',4,'FontSize',16);
%exportgraphics(fig,[test_str,'_Charge_error_coarse_PIF_',grid_str,'_',Np_str,'_dt_003125.pdf']);
exportgraphics(fig,[test_str,'_charge_error_',grid_str,'_',Np_str,'.pdf']);

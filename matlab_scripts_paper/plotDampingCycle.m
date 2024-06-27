clear all;
kx = 0.5;
Lx = 2 * pi / kx;
test = 'linear';

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
test_str = 'LandauDamping';
dir = ['../data/PinT/', test_str,'/corrected_shape_function/Conservation_studies/',time_str,'/', Np_str, '/', grid_str, '/'];
iterRank = cell(nranks, nCycles);
        
fig=figure;
for nc=1:nCycles
    for r=1:nranks
        file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/para_tol_1em8/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
        B = readmatrix(file,'NumHeaderLines',1,'Delimiter',' ');
        iterRank{r,nc} = B(:,1);
    end
end

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
    
    EzIter = cell(max_iter, 1);
    timeIter   = cell(max_iter, 1);
    for iter=1:max_iter
        shift=1;
        total = 0;
        timeIter{iter} = zeros(1,1);
        energyIter{iter} = zeros(1,3);
        for r=start_rank:step:end_rank
            r
            if(iter <= iterRank{r,nc}(end))
                file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/para_tol_1em8/data/FieldBumponTail_rank_', num2str(r-1),'_nc_',num2str(nc),'_iter_',num2str(iter),'.csv'];
                B = readmatrix(file,'NumHeaderLines',0,'Delimiter',' ');
                total = total + size(B,1);
                EzIter{iter}(shift:total, :) = B(:,2:end);
                timeIter{iter}(shift:total, :) = B(:,1);
                shift = shift + size(B,1);
            end
        end
    end



    if(strcmp(test,'nonlinear'))
        [val1,ind1] = min(abs(timeIter{2}(:)-2.5));
        [val2,ind2] = min(abs(timeIter{2}(:)-20.592));
        gamma1 = -0.562;
        gamma2 = 0.168;
        theo_rate2=exp(gamma2 * timeIter{2}(:));
        theo_rate2=(EzIter{2}(ind2,1)/theo_rate2(ind2))*theo_rate2;
        theo_amp2=exp(0.5 * gamma2 * timeIter{2}(:));
        theo_amp2=(EzIter{2}(ind2,2)/theo_amp2(ind2))*theo_amp2;
    elseif(strcmp(test,'linear'))
        [val1,ind1] = min(abs(timeIter{2}(:)-2.5));
        gamma1 = -0.3066;
    end

    if(nc == 1)
        time1  = (0:0.05:19.2)';
        [valt,indt] = min(abs(time1(:)-2.5));
        theo_rate1=exp(gamma1 * time1);
        theo_rate1=(EzIter{2}(ind1,1)/theo_rate1(indt))*theo_rate1;
    end


    if(strcmp(test,'nonlinear'))
        if(nc == 2)
            time2  = (19.2:0.05:38.4)';
            [valt,indt] = min(abs(time2(:)-20.592));
            theo_rate2=exp(gamma2 * time2);
            theo_rate2=(EzIter{2}(ind2,1)/theo_rate2(indt))*theo_rate2;
        end
    end

    
    %set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
    %                       .75 0 .75; .75 .75 0; .25 .25 .25])
    color_map = get(0, 'DefaultAxesColorOrder');
    for iter=1:max_iter
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
        pl(nc,iter) = semilogy(timeIter{iter}(:),EzIter{iter}(:,1),'LineStyle',linestyle,'Color',color_map(nco,:),'LineWidth',1.5);
        hold on;
    end
    theoE1 = semilogy(time1(:),theo_rate1(:),'b--','LineWidth',1.5);
    if(strcmp(test,'nonlinear'))
        if(nc == 2)
            hold on;
            theoE2 = semilogy(time2(:),theo_rate2(:),'k--','LineWidth',1.5);
        end
    end
    hold on;
end

%%Plot reference curve for serial time stepping
dir_serial = ['../data/serial_time/',test_str,'_conservation_studies/corrected_shape_function/'];
A_pif=readmatrix([dir_serial,'64_64_64_Pc_10/T_192/ngpus_',num2str(sranks),...
                  '/dt_003125/fine_tol_1em7/data/FieldBumponTail_',num2str(sranks),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
semilogy(A_pif(:,1),A_pif(:,2),'k--','LineWidth',2.0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




hold off;
grid on;
xlabel('time');
ylabel('$\int E_z^2 \mathrm{d}V$');
if(strcmp(test,'nonlinear'))
    ylim([1e-3 2e3]);
end
set(gca,'Fontsize',22);
legend('$k = 1$','$k = 2$','$k = 3$','$k = 4$','$k = 5$','$k = 6$','$k = 7$','analytical','serial','Location','south','NumColumns',4,'FontSize',22);
legend('boxoff');
exportgraphics(fig,[test_str,'_weak_damp_rate_',grid_str,'_',Np_str,'.pdf']);

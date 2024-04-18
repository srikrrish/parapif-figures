clear all;
set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

nCycles = 1;

EXPORT_GRAPHICS = true;
PLOT_R_ERROR = false;
PLOT_P_ERROR = false;
%PLOT_MAX_ERROR = true;

sranks = 4;
nranks = 16;
color_map = get(0, 'DefaultAxesColorOrder');
%grid_str = '32_cube_PIF_16_cube_PIC';
grid_str = '64_cube';
%Np_str =  '6.5e5'
nufft_tol = [0.001 0.0001 0.00001 0.000001 0.0000001];
%nufft_tol = [0.001 0.0001 0.00001 0.000001];

Np_str =  'Pc_10';
time_str = 'T_768';
time_end = 76.8;

% test_str = 'PenningTrap_096';
% time_end = 9.6;
%
test_str = 'LandauDamping';
%
% test_str = 'PenningTrap_38432';
% test_str = 'PenningTrap_LongTime';
% time_end = 38.4;
%
% test_str = 'PenningTrap_576';
% time_end = 57.6;
%
% test_str = 'PenningTrap_768';
% time_end = 76.8;

% addpath(genpath('~/data/seminar-thesis/PenningTrap_LongTime/'));
%nranks = [2 4 8 16 32 64];

maxRerror = zeros(2,1);
maxPerror = zeros(2,1);
figmax=figure;
ct=1;
for t=1:length(nufft_tol) 
    dir = ['../',test_str,'/',time_str,'/',Np_str,'/', grid_str,'/'];


    %shift=1;
    Rerror = cell(nranks, nCycles);
    Perror = cell(nranks, nCycles);
    iterRank = cell(nranks, nCycles);

    for r=1:nranks
        for nc = 1:nCycles
            %B = readmatrix(['./data/Np_',Np_str,'/',grid_str,'/',test_str,'/nranks_',num2str(nranks),'/CIC/L2Error/data/localError_',num2str(r-1),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            file = [dir, num2str(4), '_cycles/',num2str(sranks),'x',num2str(nranks),'/conv_test/1em12_para_tol_1em11/',num2str(nufft_tol(t),['%0.',num2str(abs(log10(nufft_tol(t)))),'f']),'/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
            B = readmatrix(file,'NumHeaderLines',1,'Delimiter',' ');
            % B = readmatrix(['./data/Different_grid_Pc/',test_str,'/NUFFT/',Np_str,'/',grid_str,'/4_cycles/data/localError_rank_',num2str(r-1),'_nc_',num2str(nc),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            Rerror{r,nc} = B(:,2);
            Perror{r,nc} = B(:,3);
            iterRank{r,nc} = B(:,1);
        end
    end


    pl = zeros(1,1);
    %% RERROR PLOT
    if PLOT_R_ERROR
        [val, idx] = max(max(cellfun(@max, iterRank)));
        fig=figure;
        for nc=1:nCycles
            if(mod(nc,2) == 0)
                lastRank = 1;
            else
                lastRank = nranks;
            end
            max_iter = iterRank{lastRank,nc}(end);
            RerrorIter = cell(max_iter, 1);
            ranksIter = cell(max_iter, 1);
            for iter=1:max_iter
                count = 1;
                RerrorIter{iter,1} = zeros(1,1);
                ranksIter{iter,1} = zeros(1,1);
                for r=1:nranks
                    if(iter <= iterRank{r,nc}(end))
                        RerrorIter{iter,1}(count) = Rerror{r,nc}(iter);
                        
                        if(mod(nc,2) == 0)
                            ranksIter{iter,1}(count) = abs((r-1) - (nranks-1));
                        else
                            ranksIter{iter,1}(count) = r-1;
                        end

                        count = count + 1;
                    end
                end
            end


            for iter=1:max_iter
                if(iter <= 7)
                    nco = iter;
                    linestyle = '-';
                    markerstyle = '*';
                elseif((iter > 7) && (iter <=14))
                    nco = mod(iter, 7);
                    if(nco == 0)
                        nco = 1;
                    end
                    linestyle = '--';
                    markerstyle = 's';
                elseif((iter > 14) && (iter <=21))
                    nco = mod(iter, 7);
                    if(nco == 0)
                        nco = 1;
                    end
                    linestyle = '-.';
                    markerstyle = 'o';
                else
                    nco = mod(iter, 7);
                    linestyle = '-';
                    markerstyle = 'd';
                end

                xaxis = ((ranksIter{iter,1}(:) + 1) * (time_end/(nranks*nCycles))) + ((nc - 1) * (time_end/nCycles));
                pl(nc,iter) = semilogy(xaxis(:), RerrorIter{iter,1}(:),'LineStyle',linestyle,'Marker',markerstyle,'Color',color_map(nco,:),'LineWidth',1.5);
                hold on;
            end
        end
        hold off;
        grid on;
        xlabel('time in s');
        ylabel('Local Rel. error in position');
        set(gca,'Fontsize',14);
        legend(pl(idx,:),'k = 1','k = 2','k = 3','k = 4','k = 5', 'k = 6', 'k = 7', 'k = 8',...
            'Location','northoutside','Orientation','horizontal','NumColumns',4,'FontSize',14);
        if EXPORT_GRAPHICS
            exportgraphics(fig,[test_str,'_RlocalError_Vs_Iter_',grid_str,'_CIC_Np_',Np_str,num2str(Pc(p)),'_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles), '.pdf']);
        end
    end


    %% PERROR PLOT
    if PLOT_P_ERROR
        [val, idx] = max(max(cellfun(@max, iterRank)));
        fig=figure;
        for nc=1:nCycles
            if(mod(nc,2) == 0)
                lastRank = 1;
            else
                lastRank = nranks;
            end
            max_iter = iterRank{lastRank,nc}(end);
            PerrorIter = cell(max_iter, 1);
            ranksIter = cell(max_iter, 1);
            for iter=1:max_iter
                count = 1;
                PerrorIter{iter,1} = zeros(1,1);
                ranksIter{iter,1} = zeros(1,1);
                for r=1:nranks
                    if(iter <= iterRank{r,nc}(end))
                        PerrorIter{iter,1}(count) = Perror{r,nc}(iter);
                        
                        if(mod(nc,2) == 0)
                            ranksIter{iter,1}(count) = abs((r-1) - (nranks-1));
                        else
                            ranksIter{iter,1}(count) = r-1;
                        end
                        
                        count = count + 1;
                    end
                end
            end

            for iter=1:max_iter
                if(iter <= 7)
                    nco = iter;
                    linestyle = '-';
                    markerstyle = '*';
                elseif((iter > 7) && (iter <=14))
                    nco = mod(iter, 7);
                    if(nco == 0)
                        nco = 1;
                    end
                    linestyle = '--';
                    markerstyle = 's';
                elseif((iter > 14) && (iter <=21))
                    nco = mod(iter, 7);
                    if(nco == 0)
                        nco = 1;
                    end
                    linestyle = '-.';
                    markerstyle = 'o';
                else
                    nco = mod(iter, 7);
                    linestyle = '-';
                    markerstyle = 'd';
                end

                xaxis = ((ranksIter{iter,1}(:) + 1) * (time_end/(nranks*nCycles))) + ((nc - 1) * (time_end/nCycles));

                pl(nc,iter) = semilogy(xaxis(:), PerrorIter{iter,1}(:),'LineStyle',linestyle,'Marker',markerstyle,'Color',color_map(nco,:),'LineWidth',1.5);
                hold on;

            end

        end
        hold off;
        grid on;
        xlabel('time in s');
        ylabel('Local Rel. error in velocity');
        set(gca,'Fontsize',14);
        legend(pl(idx,:),'k = 1','k = 2','k = 3','k = 4','k = 5', 'k = 6', 'k = 7', 'k = 8',...
            'Location','northoutside','Orientation','horizontal','NumColumns',4,'FontSize',14);
        if EXPORT_GRAPHICS
            exportgraphics(fig,[test_str,'_PlocalError_Vs_Iter_',grid_str,'_CIC_Np_',Np_str,num2str(Pc(p)),'_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles), '.pdf']);
        end
    end


    %% MAX ERROR
    %if PLOT_MAX_ERROR
        [max_iter, nc] = max(max(cellfun(@max, iterRank)));
        %fig = figure;
        %nc = 1;
        nc
        nc = 1;
        if(mod(nc,2) == 1)
            max_iter = iterRank{nranks,nc}(end);
        else
            max_iter = iterRank{1,nc}(end);
        end

        PerrorIter = cell(max_iter, 1);
        for iter=1:max_iter
            count = 1;
            PerrorIter{iter,1} = zeros(1,1);
            for r=1:nranks
                if(iter <= iterRank{r,nc}(end))
                    PerrorIter{iter,1}(count) = Perror{r,nc}(iter);
                    count = count + 1;
                end
            end
            maxPerror(iter,t) = max(PerrorIter{iter}(:));
        end
        RerrorIter = cell(max_iter, 1);
        for iter=1:max_iter
            count = 1;
            RerrorIter{iter,1} = zeros(1,1);
            for r=1:nranks
                if(iter <= iterRank{r,nc}(end))
                    RerrorIter{iter,1}(count) = Rerror{r,nc}(iter);
                    count = count + 1;
                end
            end
            maxRerror(iter,t) = max(RerrorIter{iter}(:));
        end

        %mdl = fitlm((1:max_iter)',log10(maxRerror(1:max_iter,t)));
        %mdl.Coefficients
        %%Constants for PenningTrap and TSI 10 and for LandauDamping 0.5
        Referrorit = ((10*nufft_tol(t)).^((1:max_iter)))./(factorial((1:max_iter)));
        %Referrorit = ((10*nufft_tol(t)).^((1:max_iter)));
        %Referrorit = nufft_tol(t).*(1:max_iter);
        Referrorit(:) = (maxRerror(2,t)/Referrorit(2)) .* Referrorit(:);


        %mdl = fitlm((1:max_iter)',log10(maxRerror(1:max_iter,t)));
        %keyboard;
        
        %mdl

        % for iter=1:max_iter
        %     maxRerror(iter) = max(RerrorIter{iter}(:));
        %     maxPerror(iter) = max(PerrorIter{iter}(:));
        % end
        pt(ct) = semilogy((1:max_iter)', maxRerror(1:max_iter,t),'LineStyle','-','Marker','*','Color',color_map(t,:),'LineWidth',1.5);
        hold on;
        %semilogy((1:max_iter)', maxPerror(1:max_iter,t),'LineStyle','-','Marker','s','Color',color_map(t,:),'LineWidth',1.5);
        % ylim([1e-6 1])
        %hold on;
        semilogy((1:max_iter)', Referrorit(:),'LineStyle','--','Marker','d','Color',color_map(t,:),'LineWidth',1.5);
        %plot((1:max_iter)',mdl.Fitted(:),'LineStyle','--','Marker','d','Color',color_map(t,:),'LineWidth',1.5);
        %semilogy(mdl,'LineStyle','--','Marker','d','Color',color_map(t,:),'LineWidth',1.5);
        %end
    %end
        ct = ct + 1;
end
set(gca,'Fontsize',14);
grid on;
% end
xlabel('Iterations');
ylabel('Rel. error');
%ylim([1e-8 10]);
%xlim([1 8]);
set(gca,'Fontsize',14);
legend([pt],'$\epsilon=10^{-3}$','$\epsilon=10^{-4}$','$\epsilon=10^{-5}$','$\epsilon=10^{-6}$','$\epsilon=10^{-7}$','FontSize',12,'Location','northeast','Box','off');
%title([num2str(allCycles(iCycles)),'  ', test_str]);
%if EXPORT_GRAPHICS
exportgraphics(figmax,[test_str,'_MaxlocalError_Vs_Iter_',grid_str,'_',Np_str,'_CIC_all_epsilon_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles), '.pdf']);


iterCons = [1 2 3];
fig = figure;
for it=1:length(iterCons)
    loglog(nufft_tol(:), maxRerror(iterCons(it),:)','LineStyle','-','Marker','*','Color',color_map(it,:),'LineWidth',1.5);
    hold on;
    %loglog(nufft_tol(:), maxPerror(iterCons(it),:)','LineStyle','-','Marker','s','Color',color_map(it,:),'LineWidth',1.5);
    hold on;
    Referror(:) = ((nufft_tol(:)).^(iterCons(it)));
    Referror(:) = (maxRerror(iterCons(it),1)/Referror(1)) .* Referror(:);

    loglog(nufft_tol(:), Referror(:),'LineStyle','--','Color',color_map(it,:),'LineWidth',1.5);
end
    %loglog(grid(:), Referror(:),'k--','LineWidth',1.5);
hold off;
set(gca,'Fontsize',16);
grid on;
%ylim([1e-8 1]); 
xlabel('NUFFT tolerance');
ylabel('Rel. error');
set(gca,'Fontsize',16);
%legend('Pos. error','Vel. error','$Ngrid^{-3.0}$','Location','best','FontSize',12);
    legend('$k=1$','$\epsilon^{1}$','$k=2$','$\epsilon^{2}$',...
    '$k=3$','$\epsilon^{3}$','Location','bestoutside','FontSize',12);
legend('boxoff');
%exportgraphics(fig,[test_str,'_globalError_Vs_grid_points_Pc_',num2str(Pc(p)),'_iter_',num2str(iterCons),'_CIC_nranks_',num2str(nranks),'.pdf']);
exportgraphics(fig,[test_str,'_MaxlocalError_Vs_epsilon_',grid_str,'_',Np_str,'_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles),'.pdf']);



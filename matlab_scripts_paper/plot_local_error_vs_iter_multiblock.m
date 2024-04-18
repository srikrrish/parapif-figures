clear all;
set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

% allCycles defines which block numbers you want to plot
allCycles = [16];

EXPORT_GRAPHICS = true;
PLOT_R_ERROR = true;
PLOT_P_ERROR = true;
PLOT_MAX_ERROR = true;

sranks = 4;
nranks = 16;
color_map = get(0, 'DefaultAxesColorOrder');
%grid_str = '32_cube_PIF_16_cube_PIC';
grid_str = '64_cube';
%Np_str =  '6.5e5';
Np_str =  'Pc_10';
time_str = 'T_192_dt_003125';
time_end = 19.2;

% test_str = 'PenningTrap_096';
% time_end = 9.6;
%
test_str = 'PenningTrap';
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

%dir = ['../', test_str,'/speedup_studies/',time_str,'/', Np_str, '/', grid_str, '/'];
dir = ['../', test_str,'/',time_str,'/', Np_str, '/', grid_str, '/'];

for iCycles=1:length(allCycles)
    nCycles = allCycles(iCycles);

    %shift=1;
    Rerror = cell(nranks, nCycles);
    Perror = cell(nranks, nCycles);
    iterRank = cell(nranks, nCycles);

    for r=1:nranks
        for nc = 1:nCycles
            %B = readmatrix(['./data/Np_',Np_str,'/',grid_str,'/',test_str,'/nranks_',num2str(nranks),'/CIC/L2Error/data/localError_',num2str(r-1),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            %file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
            file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIF/coarse_tol_0.001/coarse_dt_0.003125/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
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
        %ylim([3e-2 1]);
        set(gca,'Fontsize',16);
        legend(pl(idx,:),'k = 1','k = 2','k = 3','k = 4','k = 5', 'k = 6', 'k = 7', 'k = 8',...
            'Location','northoutside','Orientation','horizontal','NumColumns',4,'FontSize',16);
        if EXPORT_GRAPHICS
            exportgraphics(fig,[test_str,'_RlocalError_Vs_Iter_',grid_str,'_CIC_Np_',Np_str,'_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles), '.pdf']);
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
        %ylim([3e-2 1]);
        set(gca,'Fontsize',16);
        legend(pl(idx,:),'k = 1','k = 2','k = 3','k = 4','k = 5', 'k = 6', 'k = 7', 'k = 8',...
            'Location','northoutside','Orientation','horizontal','NumColumns',4,'FontSize',16);
        if EXPORT_GRAPHICS
            exportgraphics(fig,[test_str,'_PlocalError_Vs_Iter_',grid_str,'_CIC_Np_',Np_str,'_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles), '.pdf']);
        end
    end


    %% MAX ERROR
    if PLOT_MAX_ERROR
        [max_iter, nc] = max(max(cellfun(@max, iterRank)));
        maxRerror = zeros(max_iter,1);
        maxPerror = zeros(max_iter,1);
        %fig = figure;
        nc

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
            maxPerror(iter) = max(PerrorIter{iter}(:));
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
            maxRerror(iter) = max(RerrorIter{iter}(:));
        end
        % for iter=1:max_iter
        %     maxRerror(iter) = max(RerrorIter{iter}(:));
        %     maxPerror(iter) = max(PerrorIter{iter}(:));
        % end
        fig=figure;
        semilogy(1:max_iter', maxRerror(:),'b-*','LineWidth',1.5);
        hold on;
        plot(1:max_iter', maxPerror(:),'b--s','LineWidth',1.5);
        % ylim([1e-6 1])
        hold off;
        set(gca,'Fontsize',14);
        grid on;
        % end
        xlabel('Iterations');
        ylabel('Rel. error');
        %ylim([1e-8 1]);
        %xlim([1 8]);
        set(gca,'Fontsize',16);
        legend('Pos. error','Vel. error','Location','southwest','FontSize',16);
        %title([num2str(allCycles(iCycles)),'  ', test_str]);
        if EXPORT_GRAPHICS
            exportgraphics(fig,[test_str,'_MaxlocalError_Vs_Iter_',grid_str,'_CIC_Np_',Np_str,'_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles), '.pdf']);
        end
    end
end

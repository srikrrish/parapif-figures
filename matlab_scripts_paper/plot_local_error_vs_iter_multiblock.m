clear all;
set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

% allCycles defines which block numbers you want to plot
allCycles = [1];

EXPORT_GRAPHICS = true;
PLOT_R_ERROR = true;
PLOT_P_ERROR = true;
PLOT_MAX_ERROR = true;

sranks = 4;
nranks = 16;
color_map = get(0, 'DefaultAxesColorOrder');
grid_str = '64_cube';
Np_str =  'Pc_10';
time_str = 'T_192_dt_003125';
time_end = 19.2;
test_str = 'TSI';

dir = ['../../', test_str,'/corrected_shape_function/Conservation_studies/',time_str,'/', Np_str, '/', grid_str, '/'];

time_ax = (0:0.1:19.2)';
para_tol = 1e-8;
for iCycles=1:length(allCycles)
    nCycles = allCycles(iCycles);
    Rerror = cell(nranks, nCycles);
    Perror = cell(nranks, nCycles);
    iterRank = cell(nranks, nCycles);

    for r=1:nranks
        for nc = 1:nCycles
            file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/para_tol_1em8/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
            B = readmatrix(file,'NumHeaderLines',1,'Delimiter',' ');
            Rerror{r,nc} = B(:,2);
            Perror{r,nc} = B(:,3);
            iterRank{r,nc} = B(:,1);
        end
    end


    pl = zeros(1,1);
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
        semilogy(time_ax,para_tol.*ones(size(time_ax)),'k--','LineWidth',1.5);
        set(gca,'Fontsize',22);
        hold off;
        grid on;
        xlabel('time');
        %ylabel('Rel. error in Pos.');
        legend(pl(idx,:),'$k = 1$','$k = 2$','$k = 3$','$k = 4$','$k = 5$','$k = 6$', '$k = 7$','$k = 8$',...
            '$k = 9$','$k = 10$','$k = 11$','$k=12$','$k=13$','Location','northoutside','Orientation','horizontal','NumColumns',4,'FontSize',22);
        if EXPORT_GRAPHICS
            exportgraphics(fig,[test_str,'_RlocalError_Vs_Iter_',grid_str,'_',Np_str,'.pdf']);
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
        semilogy(time_ax,para_tol.*ones(size(time_ax)),'k--','LineWidth',1.5);
        set(gca,'Fontsize',22);
        hold off;
        grid on;
        xlabel('time');
        %ylabel('Rel. error in Vel.');
        legend(pl(idx,:),'$k = 1$','$k = 2$','$k = 3$','$k = 4$','$k = 5$','$k = 6$', '$k = 7$','$k = 8$',...
            '$k = 9$','$k = 10$','$k = 11$','$k=12$','$k=13$','Location','northoutside','Orientation','horizontal','NumColumns',4,'FontSize',22);
        if EXPORT_GRAPHICS
            exportgraphics(fig,[test_str,'_PlocalError_Vs_Iter_',grid_str,'_',Np_str,'.pdf']);
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
        fig=figure;
        semilogy(1:max_iter', maxRerror(:),'b-*','LineWidth',1.5);
        hold on;
        plot(1:max_iter', maxPerror(:),'b--s','LineWidth',1.5);
        hold off;
        set(gca,'Fontsize',22);
        grid on;
        xlabel('Iterations');
        ylabel('Rel. error');
        legend('Pos. error','Vel. error','Location','southwest','FontSize',22);
        if EXPORT_GRAPHICS
            exportgraphics(fig,[test_str,'_MaxlocalError_Vs_Iter_',grid_str,'_',Np_str,'.pdf']);
        end
    end
end

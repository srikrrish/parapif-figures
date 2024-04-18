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

sranks = [16 16 16];
nranks = 16;
color_map = get(0, 'DefaultAxesColorOrder');
%grid_str = '32_cube_PIF_16_cube_PIC';
grid_str = '_cube';
%Np_str =  '6.5e5';
grid = [8 16 32];
Np_str =  'Pc_640';
time_str = 'T_768';
time_end = 76.8;

% test_str = 'PenningTrap_096';
% time_end = 9.6;

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
for ng=1:length(grid)
    grid(ng)
    dir = ['../', test_str, '/',time_str,'/',Np_str, '/',num2str(grid(ng)),grid_str,'/'];


    %shift=1;
    Rerror = cell(nranks, nCycles);
    Perror = cell(nranks, nCycles);
    iterRank = cell(nranks, nCycles);

    for r=1:nranks
        for nc = 1:nCycles
            %B = readmatrix(['./data/Np_',Np_str,'/',grid_str,'/',test_str,'/nranks_',num2str(nranks),'/CIC/L2Error/data/localError_',num2str(r-1),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            file = [dir, num2str(4), '_cycles/',num2str(sranks(ng)),'x',num2str(nranks),'/coarse_PIC/fine_1em12/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
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
            exportgraphics(fig,[test_str,'_RlocalError_Vs_Iter_',num2str(grid(ng)),grid_str,'_CIC_',Np_str,'_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles), '.pdf']);
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
            exportgraphics(fig,[test_str,'_PlocalError_Vs_Iter_',num2str(grid(ng)),grid_str,'_CIC_Np_',Np_str,'_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles), '.pdf']);
        end
    end


    %% MAX ERROR
    %if PLOT_MAX_ERROR
        [max_iter, nc] = max(max(cellfun(@max, iterRank)));
        %maxRerror = zeros(max_iter,1);
        %maxPerror = zeros(max_iter,1);
        %fig = figure;
        %nc = 1;
        nc
        nc = 1;
        max_iter = iterRank{nranks,nc}(end);
        %max_iter = iterRank{1,nc}(end);

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
            maxPerror(iter,ng) = max(PerrorIter{iter}(:));
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
            maxRerror(iter,ng) = max(RerrorIter{iter}(:));
        end

        %Referrorit = (((0.5*4*pi)/grid(ng)).^(2.*(1:max_iter)))./(factorial((1:max_iter)-1));
        Referrorit = (((0.5*4*pi)/grid(ng)).^(2.*(1:max_iter)))./(factorial((1:max_iter)));
        %Referrorit = (((2.75*25)/grid(ng)).^(2.*(1:max_iter)))./(factorial((1:max_iter)));
        Referrorit(:) = (maxRerror(2,ng)/Referrorit(2)) .* Referrorit(:);

        %mdl = fitlm((1:max_iter)',log10(maxRerror(1:max_iter,ng)));
        %mdl.Coefficients
        % for iter=1:max_iter
        %     maxRerror(iter) = max(RerrorIter{iter}(:));
        %     maxPerror(iter) = max(PerrorIter{iter}(:));
        % end
        %fig=figure;
        % ylim([1e-6 1])
        pt(ct) = semilogy((1:max_iter)',maxRerror(1:max_iter,ng),'LineStyle','-','Marker','*','Color',color_map(ng,:),'LineWidth',1.5);
        %plot((1:max_iter)',log10(maxRerror(1:max_iter,ng)),'LineStyle','-','Marker','*','Color',color_map(ng,:),'LineWidth',1.5);
        hold on;
        %semilogy((1:max_iter)',maxPerror(1:max_iter,ng),'LineStyle','-','Marker','s','Color',color_map(ng,:),'LineWidth',1.5);
        %plot((1:max_iter)',log10(maxPerror(1:max_iter,ng)),'LineStyle','-','Marker','s','Color',color_map(ng,:),'LineWidth',1.5);
        %hold on;
        semilogy((1:max_iter)', Referrorit(:),'LineStyle','--','Marker','d','Color',color_map(ng,:),'LineWidth',1.5);
        %plot((1:max_iter)',log10(Referrorit(:)),'LineStyle','--','Marker','d','Color',color_map(ng,:),'LineWidth',1.5);
        %plot((1:max_iter)',mdl.Fitted(:),'LineStyle','--','Marker','d','Color',color_map(ng,:),'LineWidth',1.5);
        ct = ct + 1;
end
set(gca,'Fontsize',14);
grid on;
% end
xlabel('Iterations');
ylabel('Rel. error');
%ylim([1e-8 1]);
%xlim([1 8]);
set(gca,'Fontsize',14);
legend([pt],'$8^3$ modes','$16^3$ modes','$32^3$ modes','Location','southwest','FontSize',12);
%title([num2str(allCycles(iCycles)),'  ', test_str]);
%if EXPORT_GRAPHICS
exportgraphics(figmax,[test_str,'_MaxlocalError_Vs_Iter_all_grid_CIC_',Np_str,'_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles), '.pdf']);


%orderP = zeros(length(grid)-1,1); 
%orderR = zeros(length(grid)-1,1); 

for iter=1:4
    end_pt = length(grid); 
    %if(iter < 4)
    %    end_pt = 3;
    %elseif(iter == 4)
    %    end_pt = 3;
    %elseif(iter == 5)
    %    end_pt = 2;
    %end
    for g=1:end_pt
        if(g < end_pt)
            orderP(g,iter) = (log(maxPerror(iter,g+1)) - log(maxPerror(iter,g)))/(log(grid(g+1)) - log(grid(g)));
            orderR(g,iter) = (log(maxRerror(iter,g+1)) - log(maxRerror(iter,g)))/(log(grid(g+1)) - log(grid(g)));
        end
    end
end

orderP
orderR


iterCons = [1 2 3 4 5 6 7];
fig = figure;
for it=1:length(iterCons)
    end_pt = length(grid); 
    %if(iterCons(it) < 5)
    %    end_pt = 4;
    %%elseif(iterCons(it) == 4)
    %%    end_pt = 3;
    %elseif(iterCons(it) == 5)
    %    end_pt = 3;
    %end
    itn = it;
    if(it > 7)
        itn  = it-7; 
    end
    loglog(grid(1:end_pt), maxRerror(iterCons(it),1:end_pt)','LineStyle','-','Marker','*','Color',color_map(itn,:),'LineWidth',1.5);
    hold on;
    %loglog(grid(1:end_pt), maxPerror(iterCons(it),1:end_pt)','LineStyle','-','Marker','s','Color',color_map(itn,:),'LineWidth',1.5);
    %hold on;

    Referror(1:end_pt) = ((4*pi)./(grid(1:end_pt))).^(2*iterCons(it));
    %Referror(1:end_pt) = ((25)./(grid(1:end_pt))).^(2*iterCons(it));
    Referror(1:end_pt) = (maxRerror(iterCons(it),1)/Referror(1)) .* Referror(1:end_pt);

    loglog(grid(1:end_pt), Referror(1:end_pt),'LineStyle','--','Color',color_map(itn,:),'LineWidth',1.5);
end
    %loglog(grid(:), Referror(:),'k--','LineWidth',1.5);
hold off;
set(gca,'Fontsize',16);
grid on;
%ylim([1e-7 1]); 
xlabel('No. of modes per dimension');
ylabel('Rel. error');
set(gca,'Fontsize',16);
%legend('Pos. error, $k=1$','Vel. error, $k=1$','$Ng^{-1}$','Pos. error, $k=2$','Vel. error, $k=2$','$Ng^{-2}$',...
%'Pos. error, $k=3$','Vel. error, $k=3$','$Ng^{-3}$','Pos. error, $k=4$','Vel. error, $k=4$','$Ng^{-4}$','Pos. error, $k=5$','Vel. error, $k=5$','$Ng^{-5}$','Location','bestoutside','FontSize',12);
legend('$k=1$','$h^2$','$k=2$','$h^4$','$k=3$','$h^6$',...
'$k=4$','$h^8$','$k=5$','$h^{10}$','$k=6$','$h^{12}$','$k=7$','$h^{14}$',...
'Location','bestoutside','FontSize',12);
legend('boxoff');
%exportgraphics(fig,[test_str,'_globalError_Vs_grid_points_Pc_',num2str(Pc(p)),'_iter_',num2str(iterCons),'_CIC_nranks_',num2str(nranks),'.pdf']);
exportgraphics(fig,[test_str,'_MaxlocalError_Vs_grid_points_',Np_str,'_nranks_',num2str(nranks),'_ncycles_', num2str(nCycles),'.pdf']);


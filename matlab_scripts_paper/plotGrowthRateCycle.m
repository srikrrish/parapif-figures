clear all;

set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

nranks = 16;
sranks = 2;
nCycles = 1;
time_str = 'T_192_dt_05'
grid_str = '64_cube';
Np_str =  'Pc_10';
test_str = 'TSI';
dir = ['../', test_str,'/speedup_studies/',time_str,'/', Np_str, '/', grid_str, '/'];
iterRank = cell(nranks, nCycles);
        %t = tiledlayout(1,2,'TileSpacing','compact','Padding','none');
        %nexttile
        fig=figure;
        for nc=1:nCycles
            for r=1:nranks
                file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIF/coarse_tol_0.01/coarse_dt_0.05/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
                %file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
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
                        file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIF/coarse_tol_0.01/coarse_dt_0.05/data/FieldBumponTail_rank_', num2str(r-1),'_nc_',num2str(nc),'_iter_',num2str(iter),'.csv'];
                        %file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/data/FieldBumponTail_rank_', num2str(r-1),'_nc_',num2str(nc),'_iter_',num2str(iter),'.csv'];
                        B = readmatrix(file,'NumHeaderLines',0,'Delimiter',' ');
                        total = total + size(B,1);
                        EzIter{iter}(shift:total, :) = B(:,2:end);
                        timeIter{iter}(shift:total, :) = B(:,1);
                        shift = shift + size(B,1);
                    end
                end
            end
            
            
            if(strcmp(test_str,'BTI'))
                [val1,ind1] = min(abs(timeIter{2}(:)-15.0));
                gamma1 = 0.1779*2;
            elseif(strcmp(test_str,'TSI'))
                [val1,ind1] = min(abs(timeIter{2}(:)-15.0));
                gamma1 = 0.2476*2;
            end
            
            
            theo_rate1=exp(gamma1 * timeIter{2}(:));
            theo_rate1=(EzIter{2}(ind1,1)/theo_rate1(ind1))*theo_rate1;
            
            
            %pif = zeros(max_iter,1);
            set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                                   .75 0 .75; .75 .75 0; .25 .25 .25])
            color_map = get(0, 'DefaultAxesColorOrder');
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
                %xaxis = ((ranksIter{iter,1}(:) + 1) * (19.2/(nranks*nCycles))) + ((nc - 1) * (19.2/nCycles));
                pl(nc, iter) = semilogy(timeIter{iter}(:),EzIter{iter}(:,1),'LineStyle',linestyle,'Color',color_map(nco,:),'LineWidth',1.5);
                hold on;
            end
            if(nc == 1)
                theoE1 = semilogy(timeIter{2}(:),theo_rate1,'k--','LineWidth',1.5);
            end
            hold on;
        end
        hold off;
        grid on;
        %pbaspect([1 1 1])
        %xlim([0 50]);
        %xlim([0 50]);
        %ylim([1e-1 1e1]);
        xlabel('time');
        ylabel('$\int E_z^2 \mathrm{d}V$');
        %ylabel('$||E_z||_\infty$');
        set(gca,'Fontsize',18);
        
        %theo_rate1=exp(0.5 * gamma1 * timeIter{2}(:));
        %theo_rate1=(EzIter{2}(ind1,2)/theo_rate1(ind1))*theo_rate1;
        %
        %nexttile
        %for iter=1:max_iter
        %    if(iter <= 7)
        %        nc = iter;
        %        linestyle = '-';
        %        markerstyle = '*';
        %    elseif((iter > 7) && (iter <=14))
        %        nc = mod(iter, 7);
        %        if(nc == 0)
        %            nc = 1;
        %        end
        %        linestyle = '--';
        %        markerstyle = 's';
        %    elseif((iter > 14) && (iter <=21))
        %        nc = mod(iter, 7);
        %        if(nc == 0)
        %            nc = 1;
        %        end
        %        linestyle = '-.';
        %        markerstyle = 'o';
        %    else
        %        nc = mod(iter, 7);
        %        linestyle = '-';
        %        markerstyle = 'd';
        %    end
        %    pif(iter) = semilogy(timeIter{iter}(:),EzIter{iter}(:,2),'LineStyle',linestyle,'Color',color_map(nc,:),'LineWidth',1.0);
        %    hold on;
        %end
        %theoE1 = semilogy(timeIter{2}(:),theo_rate1,'r--','LineWidth',1.5);
        %hold off;
        %grid on;
        %pbaspect([1 1 1])
        %xlabel('time');
        %ylabel('$||E_z||_\infty$');
        %xlim([0 50]);
        %ylim([1e-3 1e-1]);
        %set(gca,'Fontsize',16);
        %leg = legend(pif,'k = 1','k = 2','k = 3','k = 4','k = 5','k = 6','k = 7','k = 8','k = 9','k = 10','k = 11','k = 12','k = 13','k = 14','k = 15','k = 16','Orientation','horizontal','NumColumns',4,'FontSize',14);
        %leg.Layout.Tile = 'north';
        %legend(pl(idx,:),'k = 1','k = 2','k = 3','k = 4','k = 5','k = 6','k = 7','k = 8','Location','northoutside','Orientation','horizontal','NumColumns',4,'FontSize',14);
        legend(pl(idx,:),'k = 1','k = 2','k = 3','k = 4','k = 5','k = 6','k = 7','k = 8','Location','southeast','FontSize',18);
        %exportgraphics(fig,[test_str,'_growth_rate_energy_grid_',num2str(grid(g)),'_cube_Pc_',num2str(Pc(p)),'_CIC_nranks_',num2str(nranks),'.pdf']);
        exportgraphics(fig,[test_str,'_growth_rate_',grid_str,'_',Np_str,'_',time_str,'_ncycles_',num2str(nCycles),'_nranks_',num2str(nranks),'.pdf']);

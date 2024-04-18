clear all;

colormap default;
%map = brewermap(8,'*Dark2');
%set(0, 'DefaultAxesColorOrder', map([1,2,4,5,6,7,3],:))
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
%dir = ['./data/space_time/', test_str, '/strong/',time_str,'/', Np_str, '/', grid_str, '/'];
%set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
%                                   .75 0 .75; .75 .75 0; .25 .25 .25])
dir = ['../', test_str,'/speedup_studies/',time_str,'/', Np_str, '/', grid_str, '/'];
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
                %file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
                file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIF/coarse_tol_0.01/coarse_dt_0.05/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
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
                        %file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/data/Energy_rank_', num2str(r-1),'_nc_',num2str(nc),'_iter_',num2str(iter),'.csv'];
                        file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIF/coarse_tol_0.01/coarse_dt_0.05/data/Energy_rank_', num2str(r-1),'_nc_',num2str(nc),'_iter_',num2str(iter),'.csv'];
                        B = readmatrix(file,'NumHeaderLines',0,'Delimiter',' ');
                        total = total + size(B,1);
                        chargeIter{iter}(shift:total, :) = B(:,5);
                        timeIter{iter}(shift:total, :) = B(:,1);
                        shift = shift + size(B,1);
                    end
                end
            end
            
            
            %pl = zeros(max_iter,1);
            %pl = zeros(2,1);
            for iter=1:max_iter
                charge_error_pif = abs(chargeIter{iter}(:,1) - initial_charge)./abs(initial_charge);
                %avg_error_pif(nt,ng,np) = sum(energy_error_pif(:))/length(energy_error_pif);
            
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
            
                pl(nc, iter) = semilogy(timeIter{iter}(:),charge_error_pif(:),'LineStyle',linestyle,'Color',color_map(nco,:),'LineWidth',1.5);
                %pl(nc,iter) = semilogy(timeIter{iter}(:),energyIter{iter}(:,3),'LineStyle',linestyle,'Color',color_map(nco,:),'LineWidth',1.5);
                hold on;
            end
        end
        hold off;
        grid on;
        xlabel('time');
        ylabel('Rel. Charge error');
        %ylabel('Total Energy');
        %ylabel('Potential Energy');
        set(gca,'Fontsize',18);
        %legend(pl(idx,:),'k = 1','k = 2','k = 3','k = 4','k = 5','k = 6','k = 7','k = 8','Location','northoutside','Orientation','horizontal','NumColumns',4,'FontSize',18);
        %legend(pl(idx,:),'k = 1','k = 2','k = 3','k = 4','k = 5','k = 6','k = 7','k = 8','Location','southeast','FontSize',18);
        legend(pl(idx,:),'k = 1','k = 2','k = 3','k = 4','k = 5','k = 6','k = 7','k = 8',...
                         'k = 9','k = 10','k = 11','k = 12','k = 13','k = 14','k = 15','k = 16',...
                         'Location','northoutside','Orientation','horizontal','NumColumns',4,'FontSize',14);
        %exportgraphics(fig,[test_str,'_Energy_error_grid_',num2str(grid(g)),'_cube_Pc_',num2str(Pc(p)),'_CIC_nranks_',num2str(nranks),'.pdf']);
        %exportgraphics(fig,[test_str,'_Energy_error_',grid_str,'_',Np_str,'_',time_str,'_nycles_',num2str(nCycles),'_nranks_',num2str(nranks),'.pdf']);
        exportgraphics(fig,[test_str,'_Charge_error_',grid_str,'_',Np_str,'_',time_str,'_nycles_',num2str(nCycles),'_nranks_',num2str(nranks),'.pdf']);
        %exportgraphics(fig,[test_str,'_max_iter_Energy_error_',grid_str,'_',Np_str,'_',time_str,'_nycles_',num2str(nCycles),'_nranks_',num2str(nranks),'.pdf']);

%fig=figure;
%for iter=1:max_iter
%semilogy(timeIter{iter}(:),energyIter{iter}(:,2),'LineStyle','-','Color',color_map(1,:),'LineWidth',1.5);
%grid on;
%xlabel('time');
%ylabel('Kinetic energy');
%set(gca,'Fontsize',16);
%exportgraphics(fig,['PenningTrap_KineticEnergy_PIF_',grid_str,'_Np_655360_ranks_',num2str(nranks),'.pdf']);
%
%fig=figure;
%semilogy(timeIter{max_iter}(:),energyIter{max_iter}(:,1),'LineStyle','-','Color',color_map(1,:),'LineWidth',1.5);
%grid on;
%xlabel('time');
%ylabel('Potential energy');
%set(gca,'Fontsize',16);
%exportgraphics(fig,['PenningTrap_PotentialEnergy_PIF_',grid_str,'_Np_655360_ranks_',num2str(nranks),'.pdf']);
%
%fig=figure;
%semilogy(timeIter{max_iter}(:),energyIter{max_iter}(:,3),'LineStyle','-','Color',color_map(1,:),'LineWidth',1.5);
%grid on;
%xlabel('time');
%ylabel('Total energy');
%set(gca,'Fontsize',16);
%exportgraphics(fig,['PenningTrap_TotalEnergy_PIF_',grid_str,'_Np_655360_ranks_',num2str(nranks),'.pdf']);


%avg_error_ideal = dtime.^2;
%error_ideal = (avg_error_pif(2,2,1)/avg_error_ideal(2)) .* avg_error_ideal(:);
%fig=figure;
%pif = zeros(length(grids)*length(particles),1);
%count=1;
%for np=1:length(particles)
%    for ng=1:length(grids)
%        pif(count) = loglog(dtime(:),avg_error_pif(:,ng,np),'LineWidth',1.5);
%        hold on;
%        count = count+1;
%    end
%end
%pif(count)=loglog(dtime(:),error_ideal(:),'k--*','LineWidth',1.0);
%hold off;
%grid on;
%xlabel('$\Delta t$');
%ylabel('Rel. Energy error');
%set(gca,'Fontsize',12);
%%legend([pif(1) pif(2) pif(3) pif(4) pif(5) pif(6) pif(7) pif(8) ideal],'$2^3$ modes, Np=$655,360$','$4^3$ modes, Np=$655,360$','$8^3$ modes, Np = $655,360$','$16^3$ modes, Np = $655,360$',...
%% '$32^3$ modes, Np=$655,360$','$2^3$ modes, Np = $83,886,080$','$4^3$ modes, Np = $83,886,080$','$8^3$ modes, Np = $83,886,080$','$\Delta t^2$','Location','southeast'); 
%%legend([pif(1) pif(2) pif(3) pif(4) pif(5) pif(6) ideal],'$2^3$ modes, Np=$655,360$','$4^3$ modes, Np=$655,360$','$8^3$ modes, Np = $655,360$','$16^3$ modes, Np = $655,360$',...
%% '$32^3$ modes, Np=$655,360$','$4^3$ modes, Np = $83,886,080$','$\Delta t^2$','Location','southeast'); 
%legend(pif,'$2^3$ modes, Np=$655,360$','$4^3$ modes, Np=$655,360$','$8^3$ modes, Np=$655,360$','$16^3$ modes, Np=$655,360$','$\Delta t^2$','Location','southeast'); 
%exportgraphics(fig,['Energy_convergence_pif_different_grids_particles_landau_',test,'.pdf']);

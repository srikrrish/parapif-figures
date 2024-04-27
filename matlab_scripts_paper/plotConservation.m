clear all;
kx = 0.5;
Lx = 2 * pi / kx;
nranks = 4;
test = 'PenningTrap';


set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex'); 

dtime  = [0.003125];
grids = [64];
particles = [10];
color_map = get(0, 'DefaultAxesColorOrder');
avg_error_pif = zeros(length(dtime),length(grids),length(particles));
dir = ['../../../ElectrostaticPIF/',test,'_conservation_studies/corrected_shape_function/'];
if(strcmp(test,'LandauDamping') || strcmp(test,'TSI'))
    initial_charge = -(4*pi)^3;
else
    initial_charge = -1562.5;
end
for ng=1:length(grids)
    for np=1:length(particles)
        pl = zeros(size(dtime));
        fig=figure;
        for nt=1:length(dtime)
            %A_pif=readmatrix([dir,num2str(grids(ng)),'_',num2str(grids(ng)),'_',num2str(grids(ng)),'_Pc_',num2str(particles(np)),'/T_192/ngpus_',num2str(nranks),'/dt_05/fine_tol_1em4/data/Energy_',num2str(nranks),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            A_pif=readmatrix([dir,num2str(grids(ng)),'_',num2str(grids(ng)),'_',num2str(grids(ng)),'_Pc_',num2str(particles(np)),'/T_192/ngpus_',num2str(nranks),'/dt_003125/fine_tol_1em7/data/Energy_',num2str(nranks),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            energy_error_pif = abs(A_pif(:,4) - A_pif(1,4))/abs(A_pif(1,4));
            charge_error_pif = abs(A_pif(:,5) - initial_charge)/abs(A_pif(1,5));
            mom_error_pif = abs(A_pif(:,6) - A_pif(1,6))/abs(A_pif(1,6));
        
            if(nt <= 7)
                nc = nt;
            else
                nc = mod(nt, 7);
            end
        
            semilogy(A_pif(:,1),energy_error_pif(:),'LineStyle','-','Color',color_map(nc,:),'LineWidth',1.5);
            hold on;  
            semilogy(A_pif(:,1),charge_error_pif(:),'LineStyle','-','Color',color_map(nc+1,:),'LineWidth',1.5);
            hold on;  
            semilogy(A_pif(:,1),mom_error_pif(:),'LineStyle','-','Color',color_map(nc+2,:),'LineWidth',1.5);
        
        end

        hold off;
        grid on;
        xlabel('time');
        ylabel('Rel. error');
        set(gca,'Fontsize',16);
        legend('Energy','Charge','Momentum','Location','best','Fontsize',16); 
        exportgraphics(fig,['Conservation_pif_',test,'_grid_',num2str(grids(ng)),'_cube_Pc_',num2str(particles(np)),'_dt_003125.pdf']);
    end
end

clear all;
%addpath('/data/user/muralikrishnan/ippl/build_gwendolen_release_cuda/test/particle/LandauDamping_PizDaint_benchmarks/nonlinear/strong_scaling/github_repo');
kx = 0.5;
Lx = 2 * pi / kx;
nranks = 4;
test = 'PenningTrap';


%map = brewermap(8,'*Dark2');
%set(0, 'DefaultAxesColorOrder', map([1,2,4,5,6,7,3],:))
set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex'); 

dtime  = [0.05];
grids = [64];
particles = [10];
color_map = get(0, 'DefaultAxesColorOrder');
avg_error_pif = zeros(length(dtime),length(grids),length(particles));
dir = ['../',test,'_strong_scaling_studies/'];
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
            %A_pif = readmatrix(['./PIF/LandauDampingPIF_energy_convergence/Np_',num2str(particles(np)),...
            %        '/',num2str(grids(ng)),'_cube/dt_',num2str(dtime(nt)),'/data/Energy_4.csv'],'NumHeaderLines',1,'Delimiter',' ');
            %A_pif=readmatrix(['./',test,'/',Np_str,'/',grid_str,'/dt_',num2str(dtime(nt)),'/data/Energy_',num2str(nranks),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            %A_pif=readmatrix(['./NUFFT/32_cube_655360/tol_1em2/data/Energy_',num2str(nranks),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            A_pif=readmatrix([dir,num2str(grids(ng)),'_',num2str(grids(ng)),'_',num2str(grids(ng)),'_Pc_',num2str(particles(np)),'/T_192/ngpus_',num2str(nranks),'/dt_05/fine_tol_1em4/test_wo_kokkos_custom/ngpus_2/data/Energy_',num2str(2),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            %energy_error_pif = abs(A_pif(:,4) - A_pif(1,4))/abs(A_pif(1,4));
            energy_error_pif = A_pif(:,4);
            charge_error_pif = abs(A_pif(:,5) - initial_charge)/abs(A_pif(1,5));
            %mom_error_pif = abs(A_pif(:,6) - A_pif(1,6))/abs(A_pif(1,6));
            mom_error_pif = A_pif(:,6);
            avg_error_pif(nt,ng,np) = sum(energy_error_pif(:))/length(energy_error_pif);
        
            if(nt <= 7)
                nc = nt;
            else
                nc = mod(nt, 7);
            end
        
            %semilogy(A_pif(:,1),energy_error_pif(:),'LineStyle','-','Color',color_map(nc,:),'LineWidth',1.5);
            %hold on;  
            %semilogy(A_pif(:,1),charge_error_pif(:),'LineStyle','-','Color',color_map(nc+1,:),'LineWidth',1.5);
            %hold on;  
            semilogy(A_pif(:,1),mom_error_pif(:),'LineStyle','-','Color',color_map(nc+2,:),'LineWidth',1.5);
        
        end

        hold off;
        grid on;
        xlabel('time');
        ylabel('Rel. error');
        set(gca,'Fontsize',14);
        legend('Energy','Charge','Momentum','Location','best'); 
        %legend(pl,'$\Delta t$ = 0.05','$\Delta t$ = 0.025','$\Delta t$ = 0.0125','$\Delta t$ = 0.1','$\Delta t$ = 0.05',...
        %       '$\Delta t$ = 0.025','Location','southeast'); 
        exportgraphics(fig,['Conservation_pif_',test,'_grid_',num2str(grids(ng)),'_cube_Pc_',num2str(particles(np)),'.pdf']);
    end
end

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

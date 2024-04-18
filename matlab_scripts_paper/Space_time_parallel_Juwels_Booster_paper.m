clear all;
Nsystems = 2;
dir_sys = cell(Nsystems, 1);
dir_sys{1} = '/p/project/ccstma/muralikrishnan1/ippl/build_with_kokkos_4_2_00_heffte_2_4_0_finufft_2_2_0_psmpi/alpine/ElectrostaticPIF/PenningTrap_strong_scaling_studies/128_128_128_Pc_10/T_192/ngpus_';
dir_sys{2} = '/p/project/ccstma/muralikrishnan1/ippl/build_with_kokkos_4_2_00_heffte_2_4_0_finufft_2_2_0_psmpi/alpine/ElectrostaticPIF/PenningTrap_strong_scaling_studies/128_128_128_Pc_10/T_192/ngpus_';
nx = [128];
ny = [128];
nz = [128];
Pc = [10];
nxval = length(nx);
nodes = cell(length(nx), Nsystems);

nodes{1,1} = [1 2 4 8 16 32 64 128 256];
nodes{1,2} = [1 2 4 8 16 32 64 128 256 512 1024 1536];


timeIdealOverall = cell(length(nx),1);
parallel_timers=['mainTimer...........'];

%ignore_timers=['dumpData............';
%               'particlesCreation...'];

%parallel_kernels={'mainTimer'};
npar = size(parallel_timers,1);

ind_parallel = zeros(npar,1);
%ind_ignore = zeros(size(ignore_timers,1),1);

Marker{1} = '-s';
Marker{2} = '--*';
color_map = get(0, 'DefaultAxesColorOrder');
set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex'); 

%colors_s{1} = map(1,:);%'r';
%colors_s{2} = map(2,:);%'m';
%colors_s{3} = map(4,:);%'k';
%colors_s{4} = map(5,:);%'b';
%colors_s{5} = map(6,:);%'g';

time_parallel_kernels_avg = cell(Nsystems,nxval);
time_parallel_kernels_min = cell(Nsystems,nxval);
time_parallel_kernels_max = cell(Nsystems,nxval);
        
speedup_avg = cell(Nsystems,nxval);
speedup_min = cell(Nsystems,nxval);
speedup_max = cell(Nsystems,nxval);
eff_avg = cell(Nsystems,nxval);
eff_min = cell(Nsystems,nxval);
eff_max = cell(Nsystems,nxval);

for ns=1:length(nx)
    for na=1:Nsystems

        devices = nodes{ns, na};
        Np_total = nx(ns)*ny(ns)*nz(ns)*Pc(ns);

        time_parallel_kernels_avg{na,ns} = zeros(length(devices),npar);
        time_parallel_kernels_min{na,ns} = zeros(length(devices),npar);
        time_parallel_kernels_max{na,ns} = zeros(length(devices),npar);
        speedup_avg{na,ns} = zeros(length(devices),npar);
        speedup_min{na,ns} = zeros(length(devices),npar);
        speedup_max{na,ns} = zeros(length(devices),npar);
        eff_avg{na,ns} = zeros(length(devices),npar);
        eff_min{na,ns} = zeros(length(devices),npar);
        eff_max{na,ns} = zeros(length(devices),npar);
        
        
       for i=1:length(devices)
            if(na == 1)
                fileID = fopen([dir_sys{na},num2str(devices(i)),'/dt_05/fine_tol_1em4/timing.dat']);
            else
                fileID = fopen([dir_sys{na},num2str(devices(i)),'/dt_003125/fine_tol_1em7/timing.dat']);
            end

            A=textscan(fileID,'%s %f %f %f %f','HeaderLines',6,'Delimiter',' ','MultipleDelimsAsOne',1);
            fclose(fileID);
            if(na == 1)
                time_max = A{3};
                time_min = A{4};
                time_avg = A{5};
            else
                time_max = A{3}*8;
                time_min = A{4}*8;
                time_avg = A{5}*8;
            end
            for ip=1:size(parallel_timers,1)
                ind_parallel(ip) = find(strcmp(A{1},parallel_timers(ip,:))); 
            end
            %for ig=1:size(ignore_timers,1)
            %    ind_ignore(ig) = find(strcmp(A{1},ignore_timers(ig,:))); 
            %end
            %time_avg(ind_parallel(1)) = time_avg(ind_parallel(1)) - sum(time_avg(ind_ignore));
            %time_min(ind_parallel(1)) = time_min(ind_parallel(1)) - sum(time_min(ind_ignore));
            %time_max(ind_parallel(1)) = time_max(ind_parallel(1)) - sum(time_max(ind_ignore));
            
            time_parallel_kernels_avg{na,ns}(i,:) = time_avg(ind_parallel);
            time_parallel_kernels_min{na,ns}(i,:) = time_min(ind_parallel);
            time_parallel_kernels_max{na,ns}(i,:) = time_max(ind_parallel);
           
       end
        
       for i=1:npar
                speedup_avg{na,ns}(:,i) = time_parallel_kernels_avg{na,ns}(1,i)./time_parallel_kernels_avg{na,ns}(:,i);
                eff_avg{na,ns}(:,i)     = (speedup_avg{na,ns}(:,i)./(devices(:)./devices(1))) * 100;

                speedup_min{na,ns}(:,i) = time_parallel_kernels_min{na,ns}(1,i)./time_parallel_kernels_min{na,ns}(:,i);
                eff_min{na,ns}(:,i)     = (speedup_min{na,ns}(:,i)./(devices(:)./devices(1))) * 100;

                speedup_max{na,ns}(:,i) = time_parallel_kernels_max{na,ns}(1,i)./time_parallel_kernels_max{na,ns}(:,i);
                eff_max{na,ns}(:,i)     = (speedup_max{na,ns}(:,i)./(devices(:)./devices(1))) * 100;
        end
        

        time_ideal = zeros(length(devices),1);
        time_ideal(1) = time_parallel_kernels_max{na,ns}(1,1);
        time_ideal(:) = time_ideal(1)./(devices(:)./devices(1));
        timeIdealOverall{na} = time_ideal(:);
    end
end
  
    
        fig=figure;
        for ns=1:length(nx)
            for na=1:Nsystems
                devices = nodes{ns, na};
                l(na) = plot(devices(:),time_parallel_kernels_max{na,ns}(:,1),Marker{ns},'Color',color_map(na,:),'MarkerSize',8,'LineWidth',1.5);
                l(na).MarkerFaceColor = l(na).Color;
                hold on;
                plot(devices(:),timeIdealOverall{na}(:),'k--','LineWidth',1.0);
            end
            hold on;
        end
        hold on;
        
        %128^3, Pc = 10, find_dt = 0.05, 
        %For Landau damping
        %Best coarse propagator: PIC, coarse_dt=0.2, 1 cycle
        ngt_landau = [16 32 64 128 256];
        tspacetime_landau = [43.9 44.2 23.6 14.4 11.3];  
        time_ideal = zeros(length(ngt_landau),1);
        time_ideal(:) = tspacetime_landau(1); 
        time_ideal(:) = time_ideal(1)./(ngt_landau(:)./ngt_landau(1));
        p1 = plot(ngt_landau(:),tspacetime_landau(:),'-d','Color',color_map(1,:),'MarkerSize',8,'LineWidth',1.5);
        p1.MarkerFaceColor = p1.Color;
        hold on;
        plot(ngt_landau(:),time_ideal(:),'k--','LineWidth',1.0);
        
        %For Two-stream instability
        %Best coarse propagator: PIC, coarse_dt=0.2, 1 cycle
        ngt_tsi = [16 32 64 128 256];
        tspacetime_tsi = [43.5 44.1 29 17.0 12.6];  
        p2 = plot(ngt_tsi(:),tspacetime_tsi(:),'-*','Color',color_map(1,:),'MarkerSize',8,'LineWidth',1.5);
        p2.MarkerFaceColor = p2.Color;
        hold on;
        
        %For Penning trap
        %Best coarse propagator: PIC, coarse_dt=0.05, 1 cycle
        ngt_penning = [16 32 64 128 256];
        tspacetime_penning = [41.1 31.8 23 14.7 13.6];  
        p3 = plot(ngt_penning(:),tspacetime_penning(:),'-h','Color',color_map(1,:),'MarkerSize',8,'LineWidth',1.5);
        p3.MarkerFaceColor = p3.Color;

        
        %128^3, Pc = 10, find_dt = 0.003125
        %For Landau damping
        %Best coarse propagator: PIF, coarse_tol=0.001, coarse_dt=0.05, 1 cycle
        ngt_landau = [64 128 256 512 1024 1536];
        tspacetime_landau = [712 540 283.4 154.4 93.5 74.2];  
        time_ideal = zeros(length(ngt_landau),1);
        time_ideal(:) = tspacetime_landau(1); 
        time_ideal(:) = time_ideal(1)./(ngt_landau(:)./ngt_landau(1));
        p4 = plot(ngt_landau(:),tspacetime_landau(:),'-d','Color',color_map(2,:),'MarkerSize',8,'LineWidth',1.5);
        p4.MarkerFaceColor = p4.Color;
        hold on;
        plot(ngt_landau(:),time_ideal(:),'k--','LineWidth',1.0);
        
        %For Two-stream instability
        %Best coarse propagator: PIF, coarse_tol=0.00001, coarse_dt=0.05, 1 cycle
        ngt_tsi = [64 128 256 512 1024 1536];
        tspacetime_tsi = [722 551.1 290.9 162.4 101 83];  
        p5 = plot(ngt_tsi(:),tspacetime_tsi(:),'-*','Color',color_map(2,:),'MarkerSize',8,'LineWidth',1.5);
        p5.MarkerFaceColor = p5.Color;
        hold on;
        
        %For Penning trap
        %Best coarse propagator: PIC, coarse_dt=0.003125, 16 cycles
        ngt_penning = [64 128 256 512 1024 1536];
        tspacetime_penning = [656.3 661.8 407.6 231.4 139.2 110.7];  
        p6 = plot(ngt_penning(:),tspacetime_penning(:),'-h','Color',color_map(2,:),'MarkerSize',8,'LineWidth',1.5);
        p6.MarkerFaceColor = p6.Color;


        hold off;
        xlabel('No. of GPUs');
        ylabel('time (s)');
        grid on;
        xlim([1 2048]);
        ylim([1 2e4]);
        legend([l(1) l(2) p1 p2 p3],'Space only parallel, $\Delta t_f=0.05$','Space only parallel, $\Delta t_f=0.003125$',...
                'Space time parallel, Landau damping','Space time parallel, Two-stream instability','Space time parallel, PenningTrap','Location','northeast','FontSize',12);
        legend('boxoff');
        set(gca,'FontSize',16);
        set (gca, 'XScale' , 'log' , 'YScale' , 'log' );
        exportgraphics(fig,['Space_time_parallel_128_cube_Pc_10_coarse_fine_dt_miniapps_combined.pdf']);

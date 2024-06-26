clear all;
Nsystems = 1;
dir_sys = cell(Nsystems, 1);
dir_sys{1} = '/p/project/ccstma/muralikrishnan1/ippl/build_with_kokkos_4_2_00_heffte_2_4_0_finufft_2_2_0_psmpi/alpine/ElectrostaticPIF/PenningTrap_strong_scaling_studies/corrected_shape_function';
nx = [128];
ny = [128];
nz = [128];
Pc = [10];
nxval = length(nx);
nodes = cell(length(nx), Nsystems);

nodes{1,1} = [1 2 4 8 16 32 64 128];

timeIdealOverall = cell(length(nx),1);
parallel_timers=['mainTimer...........'];

%ignore_timers=['dumpData............';
%'particlesCreation...'];
%ignore_timers=['dumpData............'];

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
        dir = [dir_sys{na},'/',num2str(nx(ns)),'_',num2str(ny(ns)),'_',num2str(nz(ns)),'_Pc_',num2str(Pc(ns)),'/T_192/ngpus_'];

        Np_total = nx(ns)*ny(ns)*nz(ns)*Pc(ns);
        %Nt = 384;

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
            fileID = fopen([dir,num2str(devices(i)),'/dt_003125/fine_tol_1em7/timing.dat']);
            %fileID = fopen([dir,num2str(devices(i)),'/dt_05/fine_tol_1em4/timing.dat']);

           A=textscan(fileID,'%s %f %f %f %f','HeaderLines',6,'Delimiter',' ','MultipleDelimsAsOne',1);
           fclose(fileID);
                time_max = A{3}*8;%./Nt;
                time_min = A{4}*8;%./Nt;
                time_avg = A{5}*8;%./Nt;
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
        timeIdealOverall{ns} = time_ideal(:);
    end
end
  
    
fig=figure;
for ns=1:length(nx)
    for na=1:Nsystems
        devices = nodes{ns, na};
        l(ns) = plot(devices(:),time_parallel_kernels_max{na,ns}(:,1),Marker{na},'Color',color_map(ns,:),'MarkerSize',8,'LineWidth',1.5);
        l(ns).MarkerFaceColor = l(ns).Color;
        hold on;
    end
    hold on;
    plot(devices(:),timeIdealOverall{ns}(:),'k--','LineWidth',1.0);
end

hold off;
xlabel('No. of GPUs');
ylabel('time (s)');
grid on;
set(gca,'FontSize',16);
set (gca, 'XScale' , 'log' , 'YScale' , 'log' );
%exportgraphics(fig,['Space_parallel_64_cube_Pc_10_Penning_dt_05_fine_tol_1em4.pdf']);
exportgraphics(fig,['Space_parallel_128_cube_Pc_10_Penning_dt_003125_fine_tol_1em7.pdf']);


clear all;
Nsystems = 1;
dir_sys = cell(Nsystems, 1);
dir_sys{1} = '../data/serial_time/PenningTrap_strong_scaling_studies/corrected_shape_function/128_128_128_Pc_10/T_192/ngpus_';
nx = [128];
ny = [128];
nz = [128];
Pc = [10];
nxval = length(nx);
nodes = cell(length(nx), Nsystems);

nodes{1,1} = [1 2 4 8 16 32 64 128 256 512 1024 1536];


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
            fileID = fopen([dir_sys{na},num2str(devices(i)),'/dt_003125/fine_tol_1em7/timing.dat']);

            A=textscan(fileID,'%s %f %f %f %f','HeaderLines',6,'Delimiter',' ','MultipleDelimsAsOne',1);
            fclose(fileID);
            %We only ran for 768 time steps whereas we need timing for 6144 time steps hence the factor 
            %8 multiplication
            time_max = A{3}*8;
            time_min = A{4}*8;
            time_avg = A{5}*8;
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
        

%linear shape function
%128^3, Pc = 10, find_dt = 0.003125
%For Landau damping
%Best coarse propagator: PIF, coarse_tol=0.001, coarse_dt=0.05, 1 cycle
ngt_landau = [64 128 256 512 1024 1536];
tspacetime_landau1 = [708.3 709.5 367 241 136.4 103];  
time_ideal = zeros(length(ngt_landau),1);
time_ideal(:) = tspacetime_landau1(1); 
time_ideal(:) = time_ideal(1)./(ngt_landau(:)./ngt_landau(1));
p1 = plot(ngt_landau(:),tspacetime_landau1(:),'-d','Color',color_map(2,:),'MarkerSize',8,'LineWidth',1.5);
p1.MarkerFaceColor = p1.Color;
hold on;
plot(ngt_landau(:),time_ideal(:),'k--','LineWidth',1.0);

%For Two-stream instability
%Best coarse propagator: PIF, coarse_tol=0.0001, coarse_dt=0.05, 1 cycle
ngt_tsi = [64 128 256 512 1024 1536];
tspacetime_tsi1 = [711.5 541.2 284.9 204.6 118.2 91.4];  
p2 = plot(ngt_tsi(:),tspacetime_tsi1(:),'-*','Color',color_map(2,:),'MarkerSize',8,'LineWidth',1.5);
p2.MarkerFaceColor = p2.Color;
hold on;

%For Penning trap
%Best coarse propagator: PIC, coarse_dt=0.0125, 16 cycles
ngt_penning = [64 128 256 512 1024 1536];
tspacetime_penning1 = [709.2 717.3 477.3 278 211.8 179.3];  
p3 = plot(ngt_penning(:),tspacetime_penning1(:),'-h','Color',color_map(2,:),'MarkerSize',8,'LineWidth',1.5);
p3.MarkerFaceColor = p3.Color;
hold on;

%Order 7 B-spline shape function
%128^3, Pc = 10, find_dt = 0.003125
%For Landau damping
%Best coarse propagator: PIF, coarse_tol=0.001, coarse_dt=0.05, 1 cycle
tspacetime_landau7 = [708 536.7 256 156.8 93.1 74];  
p4 = plot(ngt_landau(:),tspacetime_landau7(:),'-.d','Color',color_map(4,:),'MarkerSize',8,'LineWidth',1.5);
p4.MarkerFaceColor = p4.Color;
hold on;

%For Two-stream instability
%Best coarse propagator: PIF, coarse_tol=0.0001, coarse_dt=0.05, 1 cycle
tspacetime_tsi7 = [717.1 543.5 256 157.6 96 76.9];  
p5 = plot(ngt_tsi(:),tspacetime_tsi7(:),'-.*','Color',color_map(4,:),'MarkerSize',8,'LineWidth',1.5);
p5.MarkerFaceColor = p5.Color;
hold on;

%For Penning trap
%Best coarse propagator: PIC, coarse_dt=0.003125, 16 cycles, 2h mesh size
tspacetime_penning7 = [653.5 658.8 358.4 200.2 122 99.1];  
p6 = plot(ngt_penning(:),tspacetime_penning7(:),'-.h','Color',color_map(4,:),'MarkerSize',8,'LineWidth',1.5);
p6.MarkerFaceColor = p6.Color;

hold off;
set(gca,'FontSize',16);
xlabel('No. of GPUs');
ylabel('time (s)');
grid on;
xlim([1 2048]);
ylim([50 2e4]);
legend([l(1) p1 p2 p3 p4 p5 p6],'Space only parallel',...
        'Space-time parallel, Landau damping, order=1','Space-time parallel, Two-stream instability, order=1','Space-time parallel, Penning trap, order=1',...
        'Space-time parallel, Landau damping, order=7','Space-time parallel, Two-stream instability, order=7','Space-time parallel, Penning trap, order=7',...
        'Location','best','FontSize',12);
legend('boxoff');
set (gca, 'XScale' , 'log' , 'YScale' , 'log' );
exportgraphics(fig,['Space_time_parallel_128_cube_Pc_10_miniapps_combined.pdf']);

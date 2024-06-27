clear all;
addpath('./heatmaps');
test='LandauDamping';
ng='64';
pc='10';
ntime='16';
%nspace='2';
%fine_dt='05';
%coarse_dt=[0.4 0.2 0.1 0.05];
%coarse_tol = {'0.001';'0.01';'0.1';'PIC'};
%ncycles = {'1';'1';'1';'1'};
nspace='4';
fine_dt='003125';
coarse_dt=[0.2 0.05 0.0125 0.003125];
coarse_tol = {'0.00001';'0.0001';'0.001';'PIC'};
ncycles = {'1';'1';'1';'1'};
%%coarse_tol_mod = {'$10^{-5}$';'$10^{-4}$';'$10^{-3}$';'PIC'};
ncoarse_dt=length(coarse_dt);
ncoarse_tol = length(coarse_tol);
dir=['../data/PinT/',test,'/corrected_shape_function/speedup_studies/T_192_dt_',fine_dt,'/Pc_',pc,'/',ng,'_cube','/'];


parallel_timers=['mainTimer...........'];

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

time_parallel_kernels_avg = zeros(ncoarse_dt,ncoarse_tol);
time_parallel_kernels_min = zeros(ncoarse_dt,ncoarse_tol);
time_parallel_kernels_max = zeros(ncoarse_dt,ncoarse_tol);

for nt=1:ncoarse_dt
    for ntol=1:ncoarse_tol
        if(ntol < ncoarse_tol)
            direc = [dir,ncycles{ntol},'_cycles','/',nspace,'x',ntime,'/coarse_PIF/coarse_tol_',coarse_tol{ntol},'/coarse_dt_',num2str(coarse_dt(nt))];
        else
            direc = [dir,ncycles{ntol},'_cycles','/',nspace,'x',ntime,'/coarse_PIC/coarse_dt_',num2str(coarse_dt(nt))];
        end
        fileID = fopen([direc,'/timing.dat']);
        A=textscan(fileID,'%s %f %f %f %f','HeaderLines',6,'Delimiter',' ','MultipleDelimsAsOne',1);
        fclose(fileID);
        time_max = A{3};
        time_min = A{4};
        time_avg = A{5};
        for ip=1:size(parallel_timers,1)
            ind_parallel(ip) = find(strcmp(A{1},parallel_timers(ip,:))); 
        end
        %for ig=1:size(ignore_timers,1)
        %    ind_ignore(ig) = find(strcmp(A{1},ignore_timers(ig,:))); 
        %end
        %time_avg(ind_parallel(1)) = time_avg(ind_parallel(1)) - sum(time_avg(ind_ignore));
        %time_min(ind_parallel(1)) = time_min(ind_parallel(1)) - sum(time_min(ind_ignore));
        %time_max(ind_parallel(1)) = time_max(ind_parallel(1)) - sum(time_max(ind_ignore));
       
        time_parallel_kernels_avg(nt,ntol) = time_avg(ind_parallel);
        time_parallel_kernels_min(nt,ntol) = time_min(ind_parallel);
        time_parallel_kernels_max(nt,ntol) = time_max(ind_parallel);
    end
end
time_step_ratio = coarse_dt./0.003125; 
%time_step_ratio = coarse_dt./0.05; 
heatmap(time_parallel_kernels_max,coarse_tol,time_step_ratio,'%0.1f','Colormap','parula','TickAngle',45,'FontSize',28,'TextColor','m','ColorBar',1);
fig = gca;
set(gca,'Fontsize',28);
xlabel('Coarse propagator');
%ylabel('Coarse time step size');
%ylabel('$\Delta t_g/\Delta t_f$');
exportgraphics(fig,[test,'_speedup_heatmap_dt_',fine_dt,'_Pc_',pc,'_',ng,'_cube.pdf']);

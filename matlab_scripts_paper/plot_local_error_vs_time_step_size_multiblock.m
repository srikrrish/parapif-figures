clear all;
set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

nCycles = 1;


sranks = 4;
nranks = 16;
color_map = get(0, 'DefaultAxesColorOrder');
grid_str = '16_cube';
Np_str =  'Pc_640';
coarse_dt = [0.4 0.2 0.1 0.05 0.025];
time_str = 'T_192_dt_003125';
time_end = 19.2;
test_str = 'LandauDamping';

maxRerror = zeros(2,1);
maxPerror = zeros(2,1);
dir = ['../../',test_str,'/corrected_shape_function/Verification_studies/',time_str,'/',Np_str,'/', grid_str,'/'];
Rerror = cell(nranks, nCycles);
Perror = cell(nranks, nCycles);
iterRank = cell(nranks, nCycles);
figiter=figure;
ct = 1;
for t=1:length(coarse_dt) 
    coarse_dt(t)
    for r=1:nranks
        for nc = 1:nCycles
            %file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIF/coarse_tol_0.000001/coarse_dt_',num2str(coarse_dt(t)),'/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
            file = [dir, num2str(nCycles),'_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIC/coarse_dt_',num2str(coarse_dt(t)),'/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
            B = readmatrix(file,'NumHeaderLines',1,'Delimiter',' ');
            Rerror{r,nc} = B(:,2);
            Perror{r,nc} = B(:,3);
            iterRank{r,nc} = B(:,1);
        end
    end

    %[max_iter, nc] = max(max(cellfun(@max, iterRank)));
    nc = 1;
    max_iter = iterRank{nranks,nc}(end);
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
        maxPerror(iter,t) = max(PerrorIter{iter}(:));
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
        maxRerror(iter,t) = max(RerrorIter{iter}(:));
    end

    %%For LandauDamping 1 itself is a good factor
    Referrorit = ((coarse_dt(t)).^(2*(1:max_iter)))./(factorial((1:max_iter)));
    Referrorit(:) = (maxRerror(2,t)/Referrorit(2)) .* Referrorit(:);

    pt(ct) = semilogy((1:max_iter)', maxRerror(1:max_iter,t),'LineStyle','-','Marker','*','Color',color_map(t,:),'LineWidth',1.5);
    hold on;
    %semilogy((1:max_iter)', maxPerror(1:max_iter,t),'LineStyle','-','Marker','s','Color',color_map(t,:),'LineWidth',1.5);
    %hold on;
    semilogy((1:max_iter)', Referrorit(:),'LineStyle','--','Color',color_map(t,:),'LineWidth',1.5);
    ct  = ct + 1;
end
set(gca,'Fontsize',16);
grid on;
xlabel('Iterations');
ylabel('Rel. error');
legend([pt],'$\Delta t_g  = 0.4$','$\Delta t_g  = 0.2$','$\Delta t_g  = 0.1$','$\Delta t_g  = 0.05$','$\Delta t_g  = 0.025$','Location','northeast','FontSize',16);  
legend('boxoff');
exportgraphics(figiter,[test_str,'_MaxlocalError_Vs_Iter_',grid_str,'_',Np_str,'_all_coarse_dt_coarse_pic.pdf']);


iterCons = [1 2 3 4 5];
figtstep = figure;
ct=1;
for it=1:length(iterCons)
    pt(ct) = loglog(coarse_dt(:), maxRerror(iterCons(it),:)','LineStyle','-','Marker','*','Color',color_map(it,:),'LineWidth',1.5);
    hold on;
    %loglog(coarse_dt(:), maxPerror(iterCons(it),:)','LineStyle','-','Marker','s','Color',color_map(it,:),'LineWidth',1.5);
    %hold on;
    Referror = ((coarse_dt(:)).^(2*iterCons(it)))./(factorial(iterCons(it)));
    Referror(:) = (maxRerror(iterCons(it),1)/Referror(1)) .* Referror(:);
    loglog(coarse_dt(:), Referror(:),'LineStyle','--','Color',color_map(it,:),'LineWidth',1.5);
    ct = ct+1;
end
hold off;
set(gca,'Fontsize',16);
grid on;
xlabel('Coarse time step size');
ylabel('Rel. error');
xticks([0.025 0.05 0.1 0.2 0.4]);
set(gca,'Fontsize',16);
legend('boxoff');
legend([pt],'$k=1$','$k=2$','$k=3$','$k=4$','$k=5$','Location','southeast','FontSize',16);
%legend([pt],'$k=1$','$k=2$','$k=3$','Location','southeast','FontSize',16);
exportgraphics(figtstep,[test_str,'_MaxlocalError_Vs_coarse_dt_coarse_pic_',grid_str,'_',Np_str,'.pdf']);



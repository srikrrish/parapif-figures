clear all;
set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

nCycles = 1;
sranks = [16 16 16 16 16];
nranks = 16;
color_map = get(0, 'DefaultAxesColorOrder');
grid_str = '32_cube';
Pc = [5 10 20 40 80];
Np_str =  'Pc_';
time_str = 'T_192_dt_05';
test_str = 'LandauDamping';
time_end = 19.2;


maxRerror = zeros(2,1);
maxPerror = zeros(2,1);
Rerror = cell(nranks, nCycles);
Perror = cell(nranks, nCycles);
iterRank = cell(nranks, nCycles);
figiter=figure;
ct = 1;
for p=1:length(Pc) 
    Pc(p)
    dir = ['../data/PinT/', test_str, '/corrected_shape_function/Verification_studies/',time_str,'/',Np_str,num2str(Pc(p)),'/', grid_str,'/'];


    for r=1:nranks
        for nc = 1:nCycles
            file = [dir, num2str(nCycles), '_cycles/',num2str(sranks(p)),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
            B = readmatrix(file,'NumHeaderLines',1,'Delimiter',' ');
            Rerror{r,nc} = B(:,2);
            Perror{r,nc} = B(:,3);
            iterRank{r,nc} = B(:,1);
        end
    end


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
        maxPerror(iter,p) = max(PerrorIter{iter}(:));
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
        maxRerror(iter,p) = max(RerrorIter{iter}(:));
    end

    %For Landau damping 1.0 is a good factor
    %Referrorit = (1.0*Pc(p)).^(-((1:max_iter)./2))./(factorial((1:max_iter)));
    %factor_pc(p) = (maxRerror(4,p)/Referrorit(4));

    pt(ct) = semilogy((1:max_iter)', maxRerror(1:max_iter,p),'LineStyle','-','Marker','*','Color',color_map(p,:),'LineWidth',1.5);
    hold on;
    semilogy((1:max_iter)', maxPerror(1:max_iter,p),'LineStyle','-','Marker','s','Color',color_map(p,:),'LineWidth',1.5);
    hold on;
    %semilogy((1:max_iter)', Referrorit(:),'LineStyle','--','Color',color_map(p,:),'LineWidth',1.5);
    ct = ct + 1;
end
set(gca,'Fontsize',16);
grid on;
xlabel('Iterations');
ylabel('Rel. error');
legend([pt],'$P_c=5$','$P_c=10$','$P_c=20$','$P_c=40$','$P_c=80$','Location','southwest','FontSize',16);
legend('boxoff');
exportgraphics(figiter,[test_str,'_MaxlocalError_Vs_Iter_',grid_str,'_all_Pc.pdf']);


iterCons = [1 2 3 4 5 6 7];
figpc = figure;
ct = 1;
for it=1:length(iterCons)
    pt(ct) = loglog(Pc(:), maxRerror(iterCons(it),:)','LineStyle','-','Marker','*','Color',color_map(it,:),'LineWidth',1.5);
    hold on;
    loglog(Pc(:), maxPerror(iterCons(it),:)','LineStyle','-','Marker','s','Color',color_map(it,:),'LineWidth',1.5);
    hold on;
    Referror(:) = Pc(:).^(-(iterCons(it)./2));
    Referror(:) = (maxRerror(iterCons(it),1)/Referror(1)) .* Referror(:);
    loglog(Pc(:), Referror(:),'LineStyle','--','Color',color_map(it,:),'LineWidth',1.5);
    ct = ct + 1;
end
hold off;
set(gca,'Fontsize',16);
grid on;
xlabel('No. of particles per cell');
%ylabel('Rel. error');
xticks([5 10 20 40 80]);
legend('boxoff');
legend([pt],'$k=1$','$k=2$','$k=3$','$k=4$','$k=5$','$k=6$','$k=7$','Location','southwest','FontSize',16);
exportgraphics(figpc,[test_str,'_MaxlocalError_Vs_Pc_',grid_str,'_all_Iter.pdf']);



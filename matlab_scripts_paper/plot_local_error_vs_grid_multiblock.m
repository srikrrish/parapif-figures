clear all;
set(0,'defaultAxesFontName','serif');
set(0,'defaultLegendFontName','serif');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

nCycles = 1;
sranks = [16 16 16 16];
nranks = 16;
color_map = get(0, 'DefaultAxesColorOrder');
grid = [16 32 64 128];
grid_str = '_cube';
Np_str =  'Pc_10';
time_str = 'T_192_dt_05';
test_str = 'PenningTrap';
time_end = 19.2;

if(strcmp(test_str,'LandauDamping') || strcmp(test_str,'TSI'))
    L = 4*pi;
else
    L = 25;
end

mesh_sizes = L./grid(:)

maxRerror = zeros(2,1);
maxPerror = zeros(2,1);
Rerror = cell(nranks, nCycles);
Perror = cell(nranks, nCycles);
iterRank = cell(nranks, nCycles);
figiter=figure;
ct=1;
for ng=1:length(grid)
    grid(ng)
    dir = ['../../', test_str, '/corrected_shape_function/Verification_studies/',time_str,'/',Np_str, '/',num2str(grid(ng)),grid_str,'/'];


    for r=1:nranks
        for nc = 1:nCycles
            file = [dir, num2str(nCycles), '_cycles/',num2str(sranks(ng)),'x',num2str(nranks),'/coarse_PIC/coarse_dt_0.05/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
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
            maxPerror(iter,ng) = max(PerrorIter{iter}(:));
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
            maxRerror(iter,ng) = max(RerrorIter{iter}(:));
        end


        %For Landau damping 0.5 is a good factor
        if(ng==1)
            factor = 0.5;
        else
            factor = 1.5;
        end
        Referrorit = ((factor.*mesh_sizes(ng)).^(2.*(1:max_iter)))./(factorial((1:max_iter)));
        Referrorit(:) = (maxRerror(8,ng)/Referrorit(8)) .* Referrorit(:);

        pt(ct) = semilogy((1:max_iter)',maxRerror(1:max_iter,ng),'LineStyle','-','Marker','*','Color',color_map(ng,:),'LineWidth',1.5);
        hold on;
        %semilogy((1:max_iter)',maxPerror(1:max_iter,ng),'LineStyle','-','Marker','s','Color',color_map(ng,:),'LineWidth',1.5);
        %hold on;
        semilogy((1:max_iter)', Referrorit(:),'LineStyle','--','Color',color_map(ng,:),'LineWidth',1.5);
        ct = ct + 1;
end
set(gca,'Fontsize',16);
grid on;
xlabel('Iterations');
ylabel('Rel. error');
legend([pt],'$h=1.56$','$h=0.78$','$h=0.39$','$h=0.19$','Location','southwest','FontSize',16);
legend('boxoff');
exportgraphics(figiter,[test_str,'_MaxlocalError_Vs_Iter_all_grid_',Np_str,'.pdf']);



iterCons = [1 2 3 4 5 6 7];
figgrid = figure;
ct=1;
for it=1:length(iterCons)
    pt(ct) = loglog(mesh_sizes(:), maxRerror(iterCons(it),:)','LineStyle','-','Marker','*','Color',color_map(it,:),'LineWidth',1.5);
    hold on;
    %loglog(mesh_sizes(:), maxPerror(iterCons(it),:)','LineStyle','-','Marker','s','Color',color_map(it,:),'LineWidth',1.5);
    %hold on;

    Referror = ((1.5.*mesh_sizes(:)).^(2.*iterCons(it)))./(factorial(iterCons(it)));
    Referror(:) = (maxRerror(iterCons(it),1)/Referror(1)) .* Referror(:);

    loglog(mesh_sizes(:), Referror(:),'LineStyle','--','Color',color_map(it,:),'LineWidth',1.5);
    ct = ct + 1;
end
hold off;
set(gca,'Fontsize',16);
grid on;
xlabel('$h$');
ylabel('Rel. error');
%legend('Pos. error, $k=1$','Vel. error, $k=1$','$Ng^{-1}$','Pos. error, $k=2$','Vel. error, $k=2$','$Ng^{-2}$',...
%'Pos. error, $k=3$','Vel. error, $k=3$','$Ng^{-3}$','Pos. error, $k=4$','Vel. error, $k=4$','$Ng^{-4}$','Pos. error, $k=5$','Vel. error, $k=5$','$Ng^{-5}$','Location','bestoutside','FontSize',12);
%legend('$k=1$','$h^2$','$k=2$','$h^4$','$k=3$','$h^6$',...
%'$k=4$','$h^8$','$k=5$','$h^{10}$','$k=6$','$h^{12}$','$k=7$','$h^{14}$',...
%'Location','bestoutside','FontSize',12);
legend('boxoff');
legend([pt],'$k=1$','$k=2$','$k=3$','$k=4$','$k=5$','$k=6$','$k=7$','Location','southeast','FontSize',16);
exportgraphics(figgrid,[test_str,'_MaxlocalError_Vs_grid_points_',Np_str,'.pdf']);


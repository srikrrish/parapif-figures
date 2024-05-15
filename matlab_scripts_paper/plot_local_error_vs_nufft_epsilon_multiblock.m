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
nufft_tol = [0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001];
Np_str =  'Pc_640';
time_str = 'T_192_dt_05';
test_str = 'LandauDamping';
time_end = 19.2;


maxRerror = zeros(2,1);
maxPerror = zeros(2,1);
dir = ['../../',test_str,'/corrected_shape_function/Verification_studies/',time_str,'/',Np_str,'/', grid_str,'/'];
Rerror = cell(nranks, nCycles);
Perror = cell(nranks, nCycles);
iterRank = cell(nranks, nCycles);
figiter=figure;
ct=1;
for t=1:length(nufft_tol) 
    nufft_tol(t)    
    for r=1:nranks
        for nc = 1:nCycles
            file = [dir, num2str(nCycles), '_cycles/',num2str(sranks),'x',num2str(nranks),'/coarse_PIF/coarse_tol_',num2str(nufft_tol(t),['%0.',num2str(abs(log10(nufft_tol(t)))),'f']),'/coarse_dt_0.05/data/localError_rank_', num2str(r-1),'_nc_',num2str(nc),'.csv'];
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

    %%For LandauDamping 2 is a good factor
    %Referrorit = ((2*nufft_tol(t)).^((1:max_iter)))./(factorial((1:max_iter)));
    %Referrorit(:) = (maxRerror(2,t)/Referrorit(2)) .* Referrorit(:);
    
    pt(ct) = semilogy((1:max_iter)', maxRerror(1:max_iter,t),'LineStyle','-','Marker','*','Color',color_map(t,:),'LineWidth',1.5);
    hold on;
    %semilogy((1:max_iter)', maxPerror(1:max_iter,t),'LineStyle','-','Marker','s','Color',color_map(t,:),'LineWidth',1.5);
    %hold on;
    %semilogy((1:max_iter)', Referrorit(:),'LineStyle','--','Color',color_map(t,:),'LineWidth',1.5);
    ct = ct + 1;
end
set(gca,'Fontsize',16);
grid on;
xlabel('Iterations');
ylabel('Rel. error');
legend([pt],'$\varepsilon=10^{-1}$','$\varepsilon=10^{-2}$','$\varepsilon=10^{-3}$','$\varepsilon=10^{-4}$','$\varepsilon=10^{-5}$','$\varepsilon=10^{-6}$','$\varepsilon=10^{-7}$','Location','northeast','FontSize',16);
legend('boxoff');
exportgraphics(figiter,[test_str,'_MaxlocalError_Vs_Iter_',grid_str,'_',Np_str,'_all_epsilon.pdf']);


iterCons = [1 2 3];
figtol = figure;
ct=1;
for it=1:length(iterCons)
    pt(ct) = loglog(nufft_tol(:), maxRerror(iterCons(it),:)','LineStyle','-','Marker','*','Color',color_map(it,:),'LineWidth',1.5);
    hold on;
    %loglog(nufft_tol(:), maxPerror(iterCons(it),:)','LineStyle','-','Marker','s','Color',color_map(it,:),'LineWidth',1.5);
    %hold on;
    Referror = nufft_tol(:).^iterCons(it);
    Referror(:) = (maxRerror(iterCons(it),1)/Referror(1)) .* Referror(:);

    loglog(nufft_tol(:), Referror(:),'LineStyle','--','Color',color_map(it,:),'LineWidth',1.5);
    ct = ct + 1;
end
hold off;
set(gca,'Fontsize',16);
grid on;
xlabel('NUFFT tolerance');
%ylabel('Rel. error');
legend('boxoff');
legend([pt],'$k=1$','$k=2$','$k=3$','Location','southeast','FontSize',16);
exportgraphics(figtol,[test_str,'_MaxlocalError_Vs_epsilon_',grid_str,'_',Np_str,'.pdf']);

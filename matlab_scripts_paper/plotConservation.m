clear all;
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
    L = 25;
end

k = 2*pi/L;
w = sqrt(1 + 3*k^2);
vph = w/k;
%gamma = 0.5*sqrt(pi/2)*(vph)^3*(1/w)^2*exp(-(vph^2)/2);
gamma = 0.5*sqrt(pi/2)*(vph)^3*(1/w)*exp(-(vph^2)/2);


for ng=1:length(grids)
    for np=1:length(particles)
        pl = zeros(size(dtime));
        fig=figure;
        for nt=1:length(dtime)
            A_pif1=readmatrix([dir,num2str(grids(ng)),'_',num2str(grids(ng)),'_',num2str(grids(ng)),'_Pc_',num2str(particles(np)),'/T_192/ngpus_',num2str(nranks),'/dt_05/fine_tol_1em4/data/Energy_',num2str(nranks),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            A_pif2=readmatrix([dir,num2str(grids(ng)),'_',num2str(grids(ng)),'_',num2str(grids(ng)),'_Pc_',num2str(particles(np)),'/T_192/ngpus_',num2str(nranks),'/dt_003125/fine_tol_1em7/data/Energy_',num2str(nranks),'.csv'],'NumHeaderLines',1,'Delimiter',' ');
            %energy_error_pif = abs(A_pif(:,4) - A_pif(1,4))/abs(A_pif(1,4));
            %charge_error_pif = abs(A_pif(:,5) - initial_charge)/abs(A_pif(1,5));
            %mom_error_pif = abs(A_pif(:,6) - A_pif(1,6))/abs(A_pif(1,6));
        
            if(nt <= 7)
                nc = nt;
            else
                nc = mod(nt, 7);
            end
        
            %semilogy(A_pif(:,1),energy_error_pif(:),'LineStyle','-','Color',color_map(nc,:),'LineWidth',1.5);
            %hold on;  
            %semilogy(A_pif(:,1),charge_error_pif(:),'LineStyle','-','Color',color_map(nc+1,:),'LineWidth',1.5);
            %hold on;  
            %semilogy(A_pif(:,1),mom_error_pif(:),'LineStyle','-','Color',color_map(nc+2,:),'LineWidth',1.5);
        
            %semilogy(A_pif1(:,1),A_pif1(:,4),'LineStyle','-','Color',color_map(nc,:),'LineWidth',1.5);
            %hold on;
            plot(A_pif2(:,1),A_pif2(:,2),'LineStyle','-','Color',color_map(nc,:),'LineWidth',1.5);
            hold on;
            plot(A_pif2(:,1),A_pif2(:,3),'LineStyle','-','Color',color_map(nc+1,:),'LineWidth',1.5);
            hold on;
            plot(A_pif2(:,1),A_pif2(:,4),'LineStyle','-','Color',color_map(nc+2,:),'LineWidth',1.5);
            hold on;
            %semilogy(A_pif1(:,1),A_pif1(:,6),'LineStyle','-','Color',color_map(nc,:),'LineWidth',1.5);
            %hold on;
            %semilogy(A_pif2(:,1),A_pif2(:,6),'LineStyle','--','Color',color_map(nc+1,:),'LineWidth',1.5);

        end

        %PEref = 0.6.*A_pif2(1,2).*exp(-gamma.*A_pif2(:,1));
        PEref = exp(-gamma.*A_pif2(:,1));
        %PEref = (A_pif1(1,2)/PEref(1)).*PEref(:);
        PEref = 1.08.*A_pif2(1,2).*PEref(:);
        semilogy(A_pif2(:,1),PEref(:),'k--','LineWidth',1.5);


        Energy_change = max(A_pif2(:,4))/min(A_pif2(:,4)) - 1
        hold off;
        grid on;
        xlabel('time');
        %ylabel('Rel. error');
        %xticks([0 1 2 3 4 5 10 15 20]);
        ylabel('Magnitude');
        set(gca,'Fontsize',16);
        %legend('Energy','Charge','Momentum','Location','best','Fontsize',16); 
        %legend('Energy, dt=0.05','Energy, dt=0.003125','Momentum, dt=0.05','Momentum, dt=0.003125','Location','best','Fontsize',16); 
        %legend('Momentum, dt=0.05','Momentum, dt=0.003125','Location','best','Fontsize',16); 
        %legend('Energy, dt=0.05','Energy, dt=0.003125','Location','best','Fontsize',16); 
        legend('PE','KE','TE','PEref','Location','best','Fontsize',16); 
        %legend('PE','PEref','Location','best','Fontsize',16); 
        %exportgraphics(fig,['Energy_Conservation_pif_',test,'_grid_',num2str(grids(ng)),'_cube_Pc_',num2str(particles(np)),'_dt_05_003125.pdf']);
        %exportgraphics(fig,['Energy_pif_',test,'_grid_',num2str(grids(ng)),'_cube_Pc_',num2str(particles(np)),'_dt_003125.pdf']);
        %exportgraphics(fig,['Energy_pif_',test,'_grid_',num2str(grids(ng)),'_cube_Pc_',num2str(particles(np)),'_dt_05.pdf']);
    end
end

figcyc = figure;
Np = length(A_pif2(:,1));
Ype = fft(A_pif2(:,2));
Yke = fft(A_pif2(:,3));
Yte = fft(A_pif2(:,4));

P2ke = abs(Yke/Np);
P1ke = P2ke(1:floor(Np/2)+1);
P1ke(2:end-1) = 2*P1ke(2:end-1);

P2pe = abs(Ype/Np);
P1pe = P2pe(1:floor(Np/2)+1);
P1pe(2:end-1) = 2*P1pe(2:end-1);

P2te = abs(Yte/Np);
P1te = P2te(1:floor(Np/2)+1);
P1te(2:end-1) = 2*P1te(2:end-1);

Tend = 19.2;
f = (2*pi/(Tend))*(0:(Np/2));
semilogy(f,P1pe,'b-','LineWidth',1.5);
hold on;
semilogy(f,P1ke,'r--','LineWidth',1.5);
hold on;
semilogy(f,P1te,'k-.','LineWidth',1.5);
hold on;

omegaplus = 4.875*ones(size(f));
omegaplusval = linspace(1e-2,1e4,length(f));
omegaz = 1.1*ones(size(f));
omegazval = linspace(1e-2,1e4,length(f));
loglog(omegaplus(:),omegaplusval(:),'m--','LineWidth',1.5);
hold on;
loglog(omegaz(:),omegazval(:),'g--','LineWidth',1.5);
hold off;
%title("Single-Sided Amplitude Spectrum of S(t)");
xlabel('frequency');
ylabel('Magnitude');
legend('PE','KE','TE','Mod. cyclotron frequency','Axial frequency');
xlim([0 20]);
exportgraphics(figcyc,['Energy_spectrum_',test,'_grid_',num2str(grids(ng)),'_cube_Pc_',num2str(particles(np)),'_dt_003125.pdf']);
%exportgraphics(figcyc,['Energy_spectrum_',test,'_grid_',num2str(grids(ng)),'_cube_Pc_',num2str(particles(np)),'_dt_05.pdf']);

load('Nitrite_Cal_no1.mat')

figure (1)
    plot((Conc(1:8))',(curr_peak(1:8))','*','MarkerEdgeColor','g','MarkerFaceColor','g')
    hold on
%     plot((Conc(9:10))',(curr_peak(9:10))','*')
%     hold on
    plot(xfit_3,yfit_3,'-g','LineWidth', 1.4)
    hold on
%     plot((xfit_3)',c_bands_N(:,2),'--r','LineWidth', 1.2);
%     hold on
%     plot((xfit_3)',c_bands_N(:,1), '--r','LineWidth', 1.2);
%     hold on
    ylabel('Peak Current $$[\mu A]$$','interpreter','latex','FontSize',16)
    xlabel('Nitrite Concentration $$[\mu M]$$','interpreter','latex','FontSize',16)
    title('Calibration curve for Nitrite','interpreter','latex','FontSize',16)
    set(gcf,'color','w')
    grid on
    LOD(1) = 3.3*(f.ME1)/(f.b1);
    load('Calibration_results_Nitrite_Only_no2.mat')
    
    plot((Conc(1:10))',(curr_peak(1:10))','s','MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
%     plot((Conc(11:11))',(curr_peak(11:11))','*')
%     hold on
    plot(xfit_3,yfit_3,'-b','LineWidth', 1.4)
    hold on
%     plot((xfit_3)',c_bands_N(:,2),'--r','LineWidth', 1.2);
%     hold on
%     plot((xfit_3)',c_bands_N(:,1), '--r','LineWidth', 1.2);
%     hold on
%     legend('Peak Current Raw Data','Excluded in Curve fit',sprintf('Curve fit: y = %f * x + %f', p1, p2),'FontSize',15);
    ylabel('Peak Current $$[\mu A]$$','interpreter','latex','FontSize',16)
    xlabel('Nitrite Concentration $$[\mu M]$$','interpreter','latex','FontSize',16)
    title('Calibration curve for Nitrite','interpreter','latex','FontSize',16)
    set(gcf,'color','w')
    grid on
    LOD(2) = 3.3*(f.ME1)/(f.b1);
   load('Calibration_results_Nitrite_Only_no3.mat')
    
    plot((Conc(1:8))',(curr_peak(1:8))','s','MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250])
    hold on
%     plot((Conc(11:11))',(curr_peak(11:11))','*')
%     hold on
    plot(xfit_3,yfit_3,'-','color',[0.9290 0.6940 0.1250],'LineWidth', 1.4)
    hold on
%     plot((xfit_3)',c_bands_N(:,2),'--r','LineWidth', 1.2);
%     hold on
%     plot((xfit_3)',c_bands_N(:,1), '--r','LineWidth', 1.2);
%     hold on
%     legend('Peak Current Raw Data','Excluded in Curve fit',sprintf('Curve fit: y = %f * x + %f', p1, p2),'FontSize',15);
    ylabel('Peak Current $$[\mu A]$$','interpreter','latex','FontSize',16)
    xlabel('Nitrite Concentration $$[\mu M]$$','interpreter','latex','FontSize',16)
    title('Calibration curve for Nitrite','interpreter','latex','FontSize',16)
    set(gcf,'color','w')
    grid on
 LOD(3) = 3.3*(f.ME1)/(f.b1);



     load('Albumin_results.mat')
    
    plot((Conc(1:9))',(curr_peak(1:9))','o','MarkerEdgeColor','red','MarkerFaceColor','r')
    hold on
%     plot((Conc(9:11))',(curr_peak(9:11))','*')
%     hold on
    plot(xfit_3,yfit_3,'-r','LineWidth', 1.4)
    hold on
  plot(xfit_4,yfit_4,'-r','LineWidth', 1.4)
    hold on
%     plot((xfit_3)',c_bands_N(:,2),'--r','LineWidth', 1.2);
%     hold on
%     plot((xfit_3)',c_bands_N(:,1), '--r','LineWidth', 1.2);
%     hold on
%     legend('Peak Current Raw Data','Excluded in Curve fit',sprintf('Curve fit: y = %f * x + %f', p1, p2),'FontSize',15);
    ylabel('Peak Current $$[\mu A]$$','interpreter','latex','FontSize',16)
    xlabel('Nitrite Concentration $$[\mu M]$$','interpreter','latex','FontSize',16)
    title('Calibration curve for Nitrite','interpreter','latex','FontSize',16)
    set(gcf,'color','w')
    grid on
 LOD(4) = 3.3*(f2.ME1)/(f.b1);
    
h = zeros(4,1);
h(1) = plot(NaN,NaN,'color',[0 1 0],'LineWidth', 2);
h(2) = plot(NaN,NaN,'color',[0 0 1],'LineWidth', 2);
h(3) = plot(NaN,NaN,'color',[0.9290 0.6940 0.1250],'LineWidth', 2);
h(4) = plot(NaN,NaN,'color',[1 0 0],'LineWidth', 2);
leg_1 = legend(h,'Nitrite_no1','Nitrite_no2','Nitrite_no3','Nitrite_with_Albumin_v1','FontSize', 10,'interpreter','latex');
% text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',12)
    
    
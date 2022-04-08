%% Cyclic voltammetry plots

i = [4 8 12 16];

figure (1)
j = 1;
ylabel('Volts [V]')
xlabel('Current [Amps]')
set(gca,'color','w')
% for j = 1:size(i,2)
%     load(sprintf('%duM.mat',i(j)))
%     if i(j) == 12
%         load(sprintf('%duL_v2.mat',i(j)))
%     end
    plot(c0_DPVivsE(:,1),c0_DPVivsE(:,2),'LineWidth',1.8);
    hold on
% end

    legend('4uM','8uM','12uM','16uM','FontSize',18);
    ylabel('Current (uA)')
    xlabel('Voltage (V)')
set(gcf,'color','w')
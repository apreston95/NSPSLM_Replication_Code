%plotIRF.m

function plotIRF(data,vars,nimp)

[nvars,temp]=size(vars);
R=round(nvars/2);

% Plot IRFs
figure('Name','IRFs to shock');  
for n=1:nvars;
        subplot(R,2,n);             if ndims(data)==4;   %if available, plot confidence bands
                                        grpyat = [(1:nimp)' -1*data(1:nimp,nvars+n,1,2); (nimp:-1:1)' -1*data(nimp:-1:1,nvars+n,1,3)];
                                        patch(grpyat(:,1),grpyat(:,2),[0.7 0.7 0.7],'edgecolor',[0.65 0.65 0.65]); 
                                    end
                                    hold on;
                                    plot(1:nimp,-1*data(1:nimp,nvars+n,1,1),'k-','LineWidth',1.5,'color','blue');
                                    plot(1:nimp,zeros(nimp),':k');
                                    grid on
                                    title(vars(n,:),'FontSize',14,'interpreter','latex') 
                                    ylabel('Percent','FontSize',12,'interpreter','latex')
                                    xlabel('Quarters','FontSize',12,'interpreter','latex')%xlabel('years','FontSize',12)%                       
                                    hold off;
end
    
  


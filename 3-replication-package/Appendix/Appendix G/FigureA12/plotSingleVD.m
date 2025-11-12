
function plotSingleVD(data,vars,nimp)

[nvars,temp]=size(vars);
R=round(nvars/2);

%Plot figure 
figure('Name','Fraction of FEV accounted for by shock'); 
for n=1:nvars;
        subplot(R,2,n);             if ndims(data)==4;   %if available, plot confidence bands
                                        grpyat = [(1:nimp)' data(1:nimp,n,1,2); (nimp:-1:1)' data(nimp:-1:1,n,1,3)];
                                        patch(grpyat(:,1),grpyat(:,2),[0.7 0.7 0.7],'edgecolor',[0.65 0.65 0.65]);
                                    end
                                    hold on;
                                    plot(1:nimp,data(1:nimp,n,1,1),'k-'); 
                                    plot(1:nimp,.5*ones(nimp),':k')
                                    title(vars(n,:),'FontSize',14)
                                    ylabel('FEV fraction','FontSize',12)
                                    xlabel('quarters','FontSize',12) %xlabel('years','FontSize',12)%
                                    axis([0 nimp 0 1]);
                                    hold off;
end
    



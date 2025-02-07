function my_subplots_github(nA,S,v,w,TCcorr,SMcorr,rTC,rSM)
    N = size(rTC{1},1);
    axis off
    set(gca,'Units','normalized','Position',[0 0 1 1]);

    text(0.10,0.01, 'Sources','Color','b','FontSize',12)
    text(0.32,0.01, 'ICA','Color','b','FontSize',12)
    text(0.51,0.01, 'ssBSS','Color','b','FontSize',12)
    text(0.69,0.01, 'Proposed H,Z','Color','b','FontSize',12)
    text(0.89,0.01, 'Proposed D,X','Color','b','FontSize',12)

    text(0.22,0.025, ['$\bar{\gamma}$ =' num2str(round(sum(SMcorr(2,:))/12,2),'%0.2f') ],'Color','b','FontSize',12, 'Interpreter','latex')
    text(0.31,0.025, ['$\bar{\gamma}$ =' num2str(round(sum(TCcorr(2,:))/12,2),'%0.2f') ],'Color','b','FontSize',12, 'Interpreter','latex')

    text(0.415,0.025, ['$\bar{\gamma}$ =' num2str(round(sum(SMcorr(3,:))/12,2),'%0.2f') ],'Color','b','FontSize',12, 'Interpreter','latex')
    text(0.505,0.025, ['$\bar{\gamma}$ =' num2str(round(sum(TCcorr(3,:))/12,2),'%0.2f') ],'Color','b','FontSize',12, 'Interpreter','latex')

    text(0.61,0.025, ['$\bar{\gamma}$ =' num2str(round(sum(SMcorr(4,:))/12,2),'%0.2f') ],'Color','b','FontSize',12, 'Interpreter','latex')
    text(0.70,0.025, ['$\bar{\gamma}$ =' num2str(round(sum(TCcorr(4,:))/12,2),'%0.2f') ],'Color','b','FontSize',12, 'Interpreter','latex')

    text(0.805,0.025, ['$\bar{\gamma}$ =' num2str(round(sum(SMcorr(5,:))/12,2),'%0.2f') ],'Color','b','FontSize',12, 'Interpreter','latex')
    text(0.895,0.025, ['$\bar{\gamma}$ =' num2str(round(sum(TCcorr(5,:))/12,2),'%0.2f') ],'Color','b','FontSize',12, 'Interpreter','latex')


    for i =1:nA 
        for j=1:S
            ihs = 0.01;  %initial_horizontal_shift
            ivs = 0.95;  %initial_vertical_shift
            shz = 0.058; %subplot_horizontal_size
            svs = S;  %subplot_vertical_size (more the better)
            hs  = 0.85;  %horizontal shift of subplots
            nR  = S+0.5; %No. of rows
            shifter = 0.23;

            if i==1
            %%
            hax=axes();
            imagesc(flipdim(reshape(abs(zscore(rSM{i}(j,:))),v,w),1)); 
            newPos=[hs*(mod(j-1,1)+ihs+0.00),   (1-1/nR)-(1/nR)*(fix((j-1)/1)+ivs-1),   shz,   1/svs];
            set(gca,'outer',newPos), 

            hax=axes();
            plot(zscore(rTC{i}(:,j)));  axis([0 N -3 3]);
            newPos=[hs*(mod(j-1,1)+ihs+0.04),   (1-1/nR)-(1/nR)*(fix((j-1)/1)+ivs-1),   3*shz,   1/svs];
            set(gca,'outer',newPos), 

            else

            %%
            zscore_rxSM = abs(zscore(rSM{i}(j,:)));
            hax=axes();
            imagesc(flipdim(reshape(zscore_rxSM,v,w),1));  colormap('hot')
            newPos=[hs*(mod(j-1,1)+ihs+shifter*(i-1)),   (1-1/nR)-(1/nR)*(fix((j-1)/1)+ivs-1),   shz,   1/svs];
            set(gca,'outer',newPos), 
            set(gca,'XTickLabel','')
            set(gca,'YTickLabel','')
            xlabel(['\gamma',' = ',num2str(round(SMcorr(i,j),2))],'color','r')

            hax=axes();
            plot(zscore(rTC{i}(:,j))); axis([0 N -3 3]);
            newPos=[hs*(mod(j-1,1)+ihs+shifter*(i-1)+0.04),   (1-1/nR)-(1/nR)*(fix((j-1)/1)+ivs-1),   3*shz,   1/svs];
            set(gca,'outer',newPos),
            set(gca,'XTickLabel','')
            set(gca,'YTickLabel','')
            xlabel(['\gamma',' = ',num2str(round(TCcorr(i,j),2))],'color','r')
            end
        end
    end

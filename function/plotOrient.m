function plotOrient(mincorr,outfile,network,station)
% original version was written by Dr. Ba Manh Le.
%
% -- use polarhistogram instead of rose.
% Updated 2022-10-26
% Yuechu Wu
% 12131066@mail.sustech.edu.cn

if nargin == 2
    network = [];station = [];
end

outdata = load(outfile);

% merge corr and x of two normalized cross correlation
data = [outdata(:,4:5);outdata(:,6:7)];

ic = data(:,1) > mincorr;
xs = data(ic,2); % all effective x
corrs = data(ic,1); % all effective corr
meanOri_coarse = mean(xs); % average value of coarse orientation

% To avoid the situation of mean([0 360])=180
for j = 1:10
    i1 = find(xs > meanOri_coarse+180);
    xs(i1) = xs(i1) - 360;
    i2 = find(xs < meanOri_coarse-180);
    xs(i2) = xs(i2) + 360;
    meanOri_coarse = mean(xs);
end

stdOri_coarse = std(xs); % standard deviation of coarse orientation

ii=find(xs < meanOri_coarse+stdOri_coarse & xs > meanOri_coarse-stdOri_coarse);
meanOri = mean(xs(ii)); % average value of optimized orientation
stdOri  = std(xs(ii));  % standard deviation of optimized orientation

figure(1);

subplot(2,1,1)
plot(data(:,1),data(:,2),'o','MarkerFaceColor','#beb8dc','MarkerEdgeColor','none','MarkerSize',4);
hold on;
plot(corr(ii),xs(ii),'o','MarkerFaceColor','#14517c','MarkerEdgeColor','none','MarkerSize',3);

% plot([mincorr mincorr],[0 360],'r','Linewidth',1.5);
plot([mincorr mincorr],[meanOri_coarse-180 meanOri_coarse+180],'Color','#d8383a','Linewidth',1.5);
legend({'All','Used',''},'Location','northwest');
title([network,' ',station,' Orientation = ',num2str(meanOri,'%.1f'),'  std = ',num2str(stdOri,'%.1f')],'FontSize',10)

set(gca,'FontSize',10,'XTick',0:0.2:1);
set(gca,'xlim',[0 1], 'ylim',[meanOri_coarse-180 meanOri_coarse+180]);
% set(gca,'xlim',[0 1], 'ylim',[0 360],'YTick',0:30:360);


figure(2)

%%%%% original version %%%%%
[tout,rout] = rose(deg2rad(xs(ii)),36);
polar(tout,rout); hold on

[xout,yout] = pol2cart(tout,rout);
% set(gca,'nextplot','add');
% fill(xout,yout,'r'); % this method may cause bugs when there are too many elements
for ifill = 1 : length(yout)/4
    fill(xout((ifill-1)*4+1:ifill*4),yout((ifill-1)*4+1:ifill*4),[131 169 212]/255);
end
[x_ori,y_ori] = pol2cart(deg2rad(meanOri),max(rout));
quiver(0,0,x_ori,y_ori,'r','LineWidth',1.5);
view(90,-90);


%%%%% new version %%%%%
%{
nbins = 18; % specify number of bins for polar histogram
h = polarhistogram(deg2rad(xs(ii,2)),nbins);
hold on;
polarplotarrow(meanOri,max(h.Values)); % (~,~,arrowhead_angle,arrowhead_scale)
% or just plot the straight line without arrow head
% polarplot(deg2rad([meanOri meanOri]),[0 max(h.Values)],'r','LineWidth',2);
set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
%}

return
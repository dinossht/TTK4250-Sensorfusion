
load joyridedata.mat;
subplot(1,2,1)
for i =1:200
    M=Z{i}(:,:);
    scatter(M(1,:),M(2,:));hold on;
    plot(Xgt(1,:),Xgt(2,:),'k','LineWidth',2); hold on
end

xlabel('x-pos')
ylabel('y-pos')
title('Joyride')

load task2data.mat;
subplot(1,2,2)
for i =1:100
    M=Z{i}(:,:);
    scatter(M(1,:),M(2,:));hold on;
    plot(Xgt(1,:),Xgt(2,:),'k','LineWidth',2); hold on
end

xlabel('x-pos')
ylabel('y-pos')
title('Simulation')


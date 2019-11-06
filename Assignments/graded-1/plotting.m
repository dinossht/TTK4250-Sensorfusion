figure(6); clf; hold on; grid on;
plot(xest(1,:), xest(2,:));
plot(Xgt(1,1+off:K+off), Xgt(2, 1+off:K+off));
axis('equal')
title(sprintf('posRMSE = %.3f, velRMSE = %.3f, peakPosDev = %.3f, peakVelDev = %.3f',posRMSE, velRMSE, peakPosDeviation, peakVelDeviation))
legend('estimate','true')

figure(7); clf; hold on; 
subplot(2,1,1);
plot(xest(5,:)); hold on; grid on;
plot(Xgt(5,1+off:K+off)); legend('true','estimate'); xlabel('timestep'); title('Turn rate'); ylabel('rad/s');

%figure(8); clf;
subplot(2,1,2);
plot(probhat');
grid on; legend('CV','CT'); title('Model probability'); xlabel('timestep'); ylabel('Prob.');

figure(9); clf;
subplot(2,1,1); 
plot(poserr); grid on;
ylabel('position error')
subplot(2,1,2);
plot(velerr); grid on;
ylabel('velocity error')

figure(10); clf;
subplot(3,1,1);
plot(NEES); grid on; hold on;
ylabel('NEES');
ciNEES = chi2inv([0.05, 0.95], 4);
inCI = sum((NEES >= ciNEES(1)) .* (NEES <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
title(sprintf('%.2f%% inside CI', inCI))

subplot(3,1,2);
plot(NEESpos); grid on; hold on;
ylabel('NEESpos');
ciNEES = chi2inv([0.05, 0.95], 2);
inCI = sum((NEESpos >= ciNEES(1)) .* (NEESpos <= ciNEES(2)))/K * 100;
title(sprintf('%.2f%% inside CI', inCI))

subplot(3,1,3);
plot(NEESvel); grid on; hold on;
ylabel('NEESvel');
ciNEES = chi2inv([0.05, 0.95], 2);
inCI = sum((NEESvel >= ciNEES(1)) .* (NEESvel <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
title(sprintf('%.2f%% inside CI', inCI))
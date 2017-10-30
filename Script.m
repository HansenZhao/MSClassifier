Ams = Mat2MSs(A172);
Nms = Mat2MSs(Normal);
s = MSPCA(0.01,0.05);
for m = 1:1:16
    s.addMS(Nms{m}(:,1),Nms{m}(:,2),strcat('N',num2str(m)));
end
for m = 1:1:13
    s.addMS(Ams{m}(:,1),Ams{m}(:,2),strcat('A',num2str(m)));
end
[coeff,score,lat] = s.plotPCA(0,0.95,[ones(1,16),2*ones(1,13)]);
bar(coeff(:,1:2));
% clc
% L1 = find(coeff(:,1)>max((coeff(:,1))*0.3));
% L1 = find(coeff(:,1)>max((coeff(:,1))*0.2));
% L1 = find(abs(coeff(:,1))>max(abs(coeff(:,1))*0.3));
% L1 = find(abs(coeff(:,1))>max(abs(coeff(:,1))*0.25));
% L1 = find(abs(coeff(:,1))>max(abs(coeff(:,1))*0.2));
% hold on;
% for m = 1:1:11
% text(s.pks(L1(m)),coeff(L1(m),1));
% end
% for m = 1:1:11
% text(L1(m),coeff(L1(m),1),num2str(s.pks(L1(m))));
% end
% c = lines(2);
% for m = 1:1:11
% text(L1(m),coeff(L1(m),1),num2str(s.pks(L1(m))),'FontSize',6,'FontColor',c(1,:));
% end
% for m = 1:1:11
% text(L1(m),coeff(L1(m),1),num2str(s.pks(L1(m))),'FontSize',6,'Color',c(1,:));
% end
% for m = 1:1:11
% text(L1(m),coeff(L1(m),1),num2str(s.pks(L1(m))),'FontSize',8,'Color',c(1,:));
% end
% L2 = find(abs(coeff(:,2))>max(abs(coeff(:,2))*0.2));
% for m = 1:1:19
% text(L2(m),coeff(L2(m),2),num2str(s.pks(L2(m))),'FontSize',8,'Color',c(2,:));
% end
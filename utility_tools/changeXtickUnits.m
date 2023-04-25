%% xTickUnits - change xtick units of current axis
function changeXtickUnits
xtick=(0:2:8)*60*60;
xtickHr=(xtick)./(60*60);
xtickHr=chop(xtickHr,2)';
set(gca,'xtick',xtick,'xtickLabel',num2str(xtickHr));
xlabel('Time (hours)');
end


clear;

filename = dir('*.dat');

for i=1:size(filename,1)
    %tvadata{i} = tvaloader(filename{i},'STD', [5 6]);
    tvadata{i} = tvaloader(filename(i).name,'STD');
end

parpool
parfor i=1:size(filename,1)
    [theta{i},tvamodel{i}] = tvainit(tvadata{i});
    [alpha,w,C,s,v,u0,chdetgm,mu] = tvadeal(tvamodel{i},1:length(theta{i}));
    [theta{i},theta_fix{i}] = tvafixer(theta{i},[],u0,tvamodel{i});
    theta_fix{i}(u0)=0;
    tvareport(tvadata{i},tvamodel{i},theta{i},theta_fix{i});
    [theta{i},tvamodel{i},change{i}] = tvashave(theta{i},tvamodel{i},tvadata{i},theta_fix{i},[0.5 100 0.01]);
    [theta{i},tvamodel{i},change{i}] = tvashave(theta{i},tvamodel{i},tvadata{i},theta_fix{i},[0.1  1000 0.01]);
    [theta{i},tvamodel{i},change{i}] = tvashave(theta{i},tvamodel{i},tvadata{i},theta_fix{i},[0.02 1000 0.01]);
end
close

for i=1:size(filename,1)
    if(i==1) tvalpr('OutputFixedt0.txt','',tvadata{i},tvamodel{i},theta{i},theta_fix{i}); end
    tvalpr('OutputFixedt0.txt',filename(i).name,tvadata{i},tvamodel{i},theta{i},theta_fix{i});
end




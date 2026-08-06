%% run_numerical_analysis.m
clear; clc; rng(20260804,'twister');
delta=0.05; rhoValues=0.5:0.1:3; MValues=[1 2 5 10]; Nmc=5000; K=80;
R=[];
for M=MValues
 for rho=rhoValues
  threshold=(rho-1)*delta; pairSuccess=0; pairEarly=0; falseCount=0; systemSuccess=0; pairTimes=[]; systemTimes=[];
  for mc=1:Nmc
   starts=-1+2*rand(M,2); targets=-1+2*rand(M,2); tau=inf(M,1); certPos=nan(M,2);
   if threshold>=0
    for j=1:M
     for k=0:K
      truePos=starts(j,:)+(k/K)*(targets(j,:)-starts(j,:));
      measurement=truePos+sample_disk_error(delta);
      if norm(measurement-targets(j,:))<=threshold+1e-14
       tau(j)=k; certPos(j,:)=truePos; break;
      end
     end
    end
   end
   success=isfinite(tau); pairSuccess=pairSuccess+sum(success); pairEarly=pairEarly+sum(tau(success)<K);
   if any(success)
    pairTimes=[pairTimes; tau(success)/K]; %#ok<AGROW>
    d=vecnorm(certPos(success,:)-targets(success,:),2,2); falseCount=falseCount+sum(d>rho*delta+1e-12);
   end
   if all(success), systemSuccess=systemSuccess+1; systemTimes(end+1,1)=max(tau)/K; end %#ok<SAGROW>
  end
  R=[R; M rho delta rho*delta pairSuccess/(Nmc*M) systemSuccess/Nmc pairEarly/(Nmc*M) mean_nan(pairTimes) mean_nan(systemTimes) falseCount/(Nmc*M) Nmc K]; %#ok<AGROW>
 end
end
names={'M','rho','delta','delta_star','pair_trigger_probability','system_trigger_probability','early_cert_probability','mean_pair_normalized_time_conditional','mean_system_normalized_time_conditional','false_cert_probability','Nmc','K'};
T=array2table(R,'VariableNames',names); writetable(T,'numerical_results.csv');
plot_numerical_results(T,MValues);
fprintf('Generated numerical_results.csv and numerical_tradeoffs.pdf/png\n');
function y=mean_nan(x), if isempty(x), y=NaN; else, y=mean(x); end, end

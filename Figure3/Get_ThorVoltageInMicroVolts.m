function [data]=Get_ThorVoltageInMicroVolts(rawdata)
% this function converts the voltage of the LFPs in microVolts to have the
% same scale as we had on Barney. 
% yvz, Aug, 15, 2023

data=rawdata;
temp=[];

for t=1:length(data.trial);

    temp{t}=[data.trial{t}(1:192,:)*1e6; data.trial{t}(193:end,:)];
    %  % Multiply by 1e6 so that everything is in microVolts similar to Barney
    sum(sum(data.trial{t}(1:192,:)*1e6));

end

data.trial=temp;
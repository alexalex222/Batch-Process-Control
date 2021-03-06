function [] = Plot_Traj(obj, varargin)

if(mod(size(varargin,2),2) ~=0)
        error('There is a problem with the passed variables');
end

if(strcmp(varargin{1},'struct'))
    user_data = varargin{2};
else
    nvarpass = size(varargin,2)/2;
    user_data = cell(size(obj.Types,2),1);
    for n = 1:nvarpass
        type = varargin{(n-1)*2+1};
        typidx = getTypeIndex(obj, type);
        user_data{typidx} = varargin{2*n};
    end
end

max_lag = max(max(obj.Lags));
mod_traj = obj.Predict_Forward(max_lag, obj.nFE - max_lag, 'struct', user_data);

num_predicted = 0;
for typidx = 1:size(obj.Types,2)
    num_predicted = num_predicted + size(obj.Predicted{typidx},2);
end

plt_idx =1;

for typidx = 1:size(obj.Types,2)
    for i = 1: size(obj.Predicted{typidx},2)
       varidx = obj.Predicted{typidx}(i);
       
       subplot(num_predicted, 1, plt_idx);
       
       plot((max_lag+1:obj.nFE).*(obj.samt), mod_traj{typidx}(varidx,:),'r')
       hold on;
       plot((1:size(user_data{typidx}(varidx,:),2)).*(obj.samt), user_data{typidx}(varidx,:),'--k');
       
       ylabel([obj.Types{typidx} ' ' num2str(varidx)]);
       xlabel('Time');
       plt_idx = plt_idx +1;
    end
end

end

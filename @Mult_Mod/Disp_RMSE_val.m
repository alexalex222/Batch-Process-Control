function [] = Disp_RMSE_val(obj)

fprintf('-----------------------------------\n');
fprintf('RMSE for validation data:\n');
for typidx = 1:size(obj.Types,2)
   for i = 1:size(obj.Predicted{typidx},2);
       varidx = obj.Predicted{typidx}(i);
       
       fprintf(['  ' obj.Types{typidx} ' ' num2str(varidx) ': \t'...
           num2str(obj.RMSE_val{typidx}(varidx)) '\n']);
   end
end

end

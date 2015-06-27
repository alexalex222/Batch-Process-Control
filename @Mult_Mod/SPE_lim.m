function SPE_lim  = SPE_lim(obj, Prob)

if(~obj.Model_Fitted)
    error('This method can only be used after a model has been fit');
end

if(~strcmpi(obj.Regression_Type, 'PLS'))
    error('This method is only available for PLS models');
end

SPE_lim = obj.g*chi2inv(Prob, obj.h);

end
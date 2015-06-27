function Tsqrd_lim = Tsqrd_lim(obj, Prob)

if(~obj.Model_Fitted)
    error('This method can only be used after a model has been fit');
end

if(~strcmpi(obj.Regression_Type, 'PLS'))
    error('This method is only available for PLS models');
end

Tsqrd_lim = (obj.numobs-1)*(obj.numobs+1)*obj.Opt_nPC/(obj.numobs*(obj.numobs-obj.Opt_nPC))*finv(Prob,obj.Opt_nPC,obj.numobs-obj.Opt_nPC);

end
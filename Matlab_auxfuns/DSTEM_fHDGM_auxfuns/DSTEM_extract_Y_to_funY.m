function [Y_long,Y_fun,Y_mat] = DSTEM_extract_Y_to_funY(DSTEM_model,gId,yhat_flag)

if isempty(yhat_flag) == 1
    yhat_flag = 0;
end

n_st = size(DSTEM_model.stem_data.stem_varset_p.Y{1,1},1);
n_dd = size(DSTEM_model.stem_data.stem_varset_p.Y{1,1},2);
n_hh = length(DSTEM_model.stem_data.stem_varset_p.Y);
for st = 1:n_st
    for dd = 1:n_dd
        for hh = 0:(n_hh-1)
            if yhat_flag == 0
                %%% Extract observed Ys
                v(hh+1,dd) = DSTEM_model.stem_data.stem_varset_p.Y{hh+1,1}(st,dd);
            else
                %%% Extract fitted Ys
                v(hh+1,dd) = DSTEM_model.stem_EM_result.y_hat_back{1,hh+1}(st,dd);
            end
        end
    end
    y(:,st) = v(:);
end

Y_mat = y;
Y_long = y(:);
Y_split = splitapply( @(x){x}, Y_long, gId );
Y_fun = cellfun(@transpose,Y_split,'UniformOutput',false);

end
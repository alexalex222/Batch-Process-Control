function [ ] = QBMPC_CALL_EXAMPLE( batch_to_run )

%try
    % file saving stuff (select directory and file name)
%     fname_base = 'QBMPC_MMA_Shift';
    fname_base = 'QBMPC_PENICILLIN';
    folder = ['.\' fname_base '\'];
    fname = [folder fname_base '.mat'];

% glucose concentration, Biomass Conc, Penicilin Conc.
    qual_target = [0.0155 12.9304  1.3280];

    % load database of batches (so that we can start from the same initial
    % condition)
    load Database\bdb_v1_3;

    
    % load the quality model
    load Q_model_v1;
    
    % Align u 
    batches = bdb.data;
    for i=1:bdb.nb
        batches.u(:,bdb.nFE+1,i) = batches.u(:,bdb.nFE,i);
    end
    batches.unom = [batches.unom batches.unom(:,end)];
    
%     I believe we don't have the multimodel yet
    % load the multi-model
    load Mult_Mod\mod;

    if(batch_to_run>0 && batch_to_run<=bdb.nb)


        % MPC CALL
        [quality_result_out percent_imp_out x_out y_out q_out u_out u_all_out fig_handle JMPC_out time_out]  ...
            = QBMPC_PMMAv2(bdb, batches, batch_to_run, Q_model, opt_mod, qual_target);

        format = '% 9.4g';

        fprintf('--------------------------------------\n');
        fprintf('                \tConv.  \t\tnMW \t\twMW\n');
        fprintf(['Target Quality: \t' num2str(qual_target(1), format ) '\t\t' num2str(qual_target(2), format) '\t\t' num2str(qual_target(3), format) '\n']);
        fprintf(['Result Quality: \t' num2str(quality_result_out(1), format) '\t\t' num2str(quality_result_out(2), format)  '\t\t' num2str(quality_result_out(3), format) '\n']);
        fprintf(['PI Quality:     \t' num2str(batches.q(1,end, batch_to_run), format) '\t\t' num2str(batches.q(2,end, batch_to_run), format) '\t\t' num2str(batches.q(3,end, batch_to_run), format) '\n']);
        fprintf(['improvement:    \t' num2str(percent_imp_out(1), format) '%%\t\t' num2str(percent_imp_out(2), format) '%%\t\t' num2str(percent_imp_out(3), format) '%%\n']);
        fprintf('\n');




        %%%% Load results from past data to append with this batch %%%%
        if(exist(folder,'dir') == 0)
           mkdir(folder); 
        end

        if(exist(fname,'file') == 0)
            i = 1;
        else
            load(fname);
            i = i+1;
        end

        %%%% Append data from this batch to that of previous batches%%%%
        batches_run(i) = batch_to_run;

        quality_pi(i,:) = batches.q(:,end,batch_to_run)';
        quality_QBMPC(i,:) = quality_result_out;
        percent_imp(i,:) = percent_imp_out;
        JMPC(:,i) = JMPC_out(:)';

        x(:,:,i) = x_out;
        y(:,:,i) = y_out;
        u(:,:,i) = u_out;
        q(:,:,i) = q_out;
        
        time(:,i) = time_out';

        u_all{i} = u_all_out;

        % Save data
        save(fname, 'batches_run', 'quality_pi', 'quality_QBMPC', 'percent_imp', 'x', 'y', 'u', 'q', 'u_all', 'i', 'JMPC', 'time');

        %Save figure
        fig_fname = [folder fname_base '_' num2str(batch_to_run)];
        saveas(fig_handle, fig_fname);
    end
%catch exception
 %  message = ['An error has occured in ' exception.stack(1).name ' at line '...
 %      num2str(exception.stack.line) '. ' exception.message];
   
 %  fprintf(message)
   %send_text(message, '5856131397');

%end
end


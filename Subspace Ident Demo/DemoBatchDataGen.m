clear vairables;
close all;

% Set batch system options
opts = BatchSysOpts('nx', 3, 'nu', 3, 'ny', 2, 'nb', 2, 'inputVar', 'PRBS');

% % Generate the system and data
% bdb1 = BatchDataGen(opts, 1, false);

% Generate a different system with the same options and data
bdb = BatchDataGen(opts, 200, true);


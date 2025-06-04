function Version_output(filename)

%% Record the matlab, mtex versions and the date
% Save them to a file

Matlab_version=version
Mtex_version=getMTEXpref().version
time_now = datetime('now','TimeZone','local','Format','dd-MMM-y HH:mm:ss')

fid = fopen(filename,'wt');
fprintf(fid, 'Matlab version: %s\nMtex version: %s\nDate: %s',...
    Matlab_version, Mtex_version, time_now);
fclose(fid);

end
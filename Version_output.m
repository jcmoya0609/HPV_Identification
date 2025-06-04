
Matlab_version=version
Mtex_version=getMTEXpref().version
time_now = datetime('now','TimeZone','local','Format','dd-MMM-y HH:mm:ss')

fid = fopen('Version_Flag.txt','wt');
fprintf(fid, 'Matlab version: %s\nMtex version: %s\nDate: %s',...
    Matlab_version, Mtex_version, time_now);
fclose(fid);
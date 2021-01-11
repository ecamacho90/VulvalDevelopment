%copyfile Parallel_function_AbsDist_Modelv9_v1_01.m Parallel_function_AbsDist_Modelv9_v1_001.m
versionmodel = '10';
versionsubmodel = '1';

for numberfunc = 1:9

fid = fopen('Parallel_function_AbsDist_Modelv10_v1_Template.m','r')
% Read txt into cell A
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);
% Change cell A
A{1} = '% Parallel_call_function_20000';
A{5} = ['function',' ','Parallel_function_AbsDist_Modelv',versionmodel,'_v',versionsubmodel,'_00',num2str(numberfunc),'()'];
% Write cell A into txt
fid = fopen(['Parallel_function_AbsDist_Modelv',versionmodel,'_v',versionsubmodel,'_00',num2str(numberfunc),'.m'], 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end

fclose(fid);
disp('Finished')
pause(5)
end

%%

for numberfunc = 10:99

fid = fopen('Parallel_function_AbsDist_Modelv10_v1_Template.m','r')
% Read txt into cell A
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);
% Change cell A
A{1} = '% Parallel_call_function_20000';
A{5} = ['function',' ','Parallel_function_AbsDist_Modelv',versionmodel,'_v',versionsubmodel,'_0',num2str(numberfunc),'()'];
% Write cell A into txt
fid = fopen(['Parallel_function_AbsDist_Modelv',versionmodel,'_v',versionsubmodel,'_0',num2str(numberfunc),'.m'], 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end

fclose(fid);
disp('Finished')
pause(5)
end



for numberfunc = 100:500

fid = fopen('Parallel_function_AbsDist_Modelv10_v1_Template.m','r')
% Read txt into cell A
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);
% Change cell A
A{1} = '% Parallel_call_function_20000';
A{5} = ['function',' ','Parallel_function_AbsDist_Modelv',versionmodel,'_v',versionsubmodel,'_',num2str(numberfunc),'()'];
% Write cell A into txt
fid = fopen(['Parallel_function_AbsDist_Modelv',versionmodel,'_v',versionsubmodel,'_',num2str(numberfunc),'.m'], 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end

fclose(fid);
disp('Finished')
pause(5)
end

disp('All finished')
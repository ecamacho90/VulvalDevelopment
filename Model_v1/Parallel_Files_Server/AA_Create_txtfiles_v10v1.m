versionmodel = '10';
versionsubmodel = '1';
datefunc = '22_04_20';

pathtofolder = '/home/';
emailuser = 'email@direc.com';
%%
for i = 1:9
    fileID = fopen(['ABCv',versionmodel,'v',versionsubmodel,'_00',num2str(i),'.txt'],'w');
    fprintf(fileID,'#!/bin/bash \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# specify the account and partition to schedule jobs \n');
    fprintf(fileID,'#SBATCH --account=commons \n');
    fprintf(fileID,'#SBATCH --partition=commons \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# name the job \n');
    fprintf(fileID,['#SBATCH --job-name=ABC','_00',num2str(i),'_v',versionmodel,'_v',versionsubmodel,'_',datefunc,' \n']);
    fprintf(fileID,'\n');
    fprintf(fileID,'# redirect job?s standard output \n');
    fprintf(fileID,['#SBATCH --output ',pathtofolder,'ABCstdoutv',versionmodel,'v',versionsubmodel,'_00',num2str(i),'.txt \n']);
    fprintf(fileID,'\n');
    fprintf(fileID,'# redirect job?s standard error \n');
    fprintf(fileID,['#SBATCH --error ',pathtofolder,'ABCstderrv',versionmodel,'v',versionsubmodel,'_00',num2str(i),'.txt \n']);
    fprintf(fileID,'\n');
    fprintf(fileID,'# maximum amount of real time that a job can be in the running state \n');
    fprintf(fileID,'#SBATCH --time=0-02:00:00 \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# number of nodes requested for the job \n');
    fprintf(fileID,'#SBATCH --nodes=1 \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# number of CPU cores requested on a single node. Because matlab can only have at most 8 jobs in parallel, request 9 cores \n');
    fprintf(fileID,'#SBATCH --cpus-per-task=11 \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# maximum amount of physical memory that the job can use \n');
    fprintf(fileID,'#SBATCH --mem-per-cpu=2G \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# email at the beginning and end of the job \n');
    fprintf(fileID,'#SBATCH --mail-type=TIME_LIMIT_80,FAIL \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# email address to receive the notifications \n');
    fprintf(fileID,['#SBATCH --mail-user=',emailuser,' \n']);
    fprintf(fileID,'# load matlab \n');
    fprintf(fileID,'module load MATLAB/2015a \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# \n');
    fprintf(fileID,'echo "I ran on:" \n');
    fprintf(fileID,'cd $SLURM_SUBMIT_DIR \n');
    fprintf(fileID,'echo $SLURM_NODELIST \n');
    fprintf(fileID,'# \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# submit matlab job \n');
    fprintf(fileID,['srun -n 1 -c 11 matlab -nodisplay -nosplash -nodesktop -r "Parallel_function_AbsDist_Modelv',versionmodel,'_v',versionsubmodel,'_00',num2str(i),'"']);
    
    fclose(fileID);
    
    
    
end

%%

for i = 10:99
    fileID = fopen(['ABCv',versionmodel,'v',versionsubmodel,'_0',num2str(i),'.txt'],'w');
    fprintf(fileID,'#!/bin/bash \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# specify the account and partition to schedule jobs \n');
    fprintf(fileID,'#SBATCH --account=commons \n');
    fprintf(fileID,'#SBATCH --partition=commons \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# name the job \n');
    fprintf(fileID,['#SBATCH --job-name=ABC','_0',num2str(i),'_v',versionmodel,'_v',versionsubmodel,'_',datefunc,' \n']);
    fprintf(fileID,'\n');
    fprintf(fileID,'# redirect job?s standard output \n');
    fprintf(fileID,['#SBATCH --output /home/ec47/ABCstdoutv',versionmodel,'v',versionsubmodel,'_0',num2str(i),'.txt \n']);
    fprintf(fileID,'\n');
    fprintf(fileID,'# redirect job?s standard error \n');
    fprintf(fileID,['#SBATCH --error /home/ec47/ABCstderrv',versionmodel,'v',versionsubmodel,'_0',num2str(i),'.txt \n']);
    fprintf(fileID,'\n');
    fprintf(fileID,'# maximum amount of real time that a job can be in the running state \n');
    fprintf(fileID,'#SBATCH --time=0-02:00:00 \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# number of nodes requested for the job \n');
    fprintf(fileID,'#SBATCH --nodes=1 \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# number of CPU cores requested on a single node. Because matlab can only have at most 8 jobs in parallel, request 9 cores \n');
    fprintf(fileID,'#SBATCH --cpus-per-task=11 \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# maximum amount of physical memory that the job can use \n');
    fprintf(fileID,'#SBATCH --mem-per-cpu=2G \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# email at the beginning and end of the job \n');
    fprintf(fileID,'#SBATCH --mail-type=TIME_LIMIT_80,FAIL \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# email address to receive the notifications \n');
    fprintf(fileID,'#SBATCH --mail-user=ec47@rice.edu \n');
    fprintf(fileID,'# load matlab \n');
    fprintf(fileID,'module load MATLAB/2015a \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# \n');
    fprintf(fileID,'echo "I ran on:" \n');
    fprintf(fileID,'cd $SLURM_SUBMIT_DIR \n');
    fprintf(fileID,'echo $SLURM_NODELIST \n');
    fprintf(fileID,'# \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# submit matlab job \n');
    fprintf(fileID,['srun -n 1 -c 11 matlab -nodisplay -nosplash -nodesktop -r "Parallel_function_AbsDist_Modelv',versionmodel,'_v',versionsubmodel,'_0',num2str(i),'"']);
    
    fclose(fileID);
    
end

%%

for i = 100:500
    
    fileID = fopen(['ABCv',versionmodel,'v',versionsubmodel,'_',num2str(i),'.txt'],'w');
    fprintf(fileID,'#!/bin/bash \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# specify the account and partition to schedule jobs \n');
    fprintf(fileID,'#SBATCH --account=commons \n');
    fprintf(fileID,'#SBATCH --partition=commons \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# name the job \n');
    fprintf(fileID,['#SBATCH --job-name=ABC','_',num2str(i),'_v',versionmodel,'_v',versionsubmodel,'_',datefunc,' \n']);
    fprintf(fileID,'\n');
    fprintf(fileID,'# redirect job?s standard output \n');
    fprintf(fileID,['#SBATCH --output /home/ec47/ABCstdoutv',versionmodel,'v',versionsubmodel,'_',num2str(i),'.txt \n']);
    fprintf(fileID,'\n');
    fprintf(fileID,'# redirect job?s standard error \n');
    fprintf(fileID,['#SBATCH --error /home/ec47/ABCstderrv',versionmodel,'v',versionsubmodel,'_',num2str(i),'.txt \n']);
    fprintf(fileID,'\n');
    fprintf(fileID,'# maximum amount of real time that a job can be in the running state \n');
    fprintf(fileID,'#SBATCH --time=0-02:00:00 \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# number of nodes requested for the job \n');
    fprintf(fileID,'#SBATCH --nodes=1 \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# number of CPU cores requested on a single node. Because matlab can only have at most 8 jobs in parallel, request 9 cores \n');
    fprintf(fileID,'#SBATCH --cpus-per-task=11 \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# maximum amount of physical memory that the job can use \n');
    fprintf(fileID,'#SBATCH --mem-per-cpu=2G \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# email at the beginning and end of the job \n');
    fprintf(fileID,'#SBATCH --mail-type=TIME_LIMIT_80,FAIL \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# email address to receive the notifications \n');
    fprintf(fileID,'#SBATCH --mail-user=ec47@rice.edu \n');
    fprintf(fileID,'# load matlab \n');
    fprintf(fileID,'module load MATLAB/2015a \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# \n');
    fprintf(fileID,'echo "I ran on:" \n');
    fprintf(fileID,'cd $SLURM_SUBMIT_DIR \n');
    fprintf(fileID,'echo $SLURM_NODELIST \n');
    fprintf(fileID,'# \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'# submit matlab job \n');
    fprintf(fileID,['srun -n 1 -c 11 matlab -nodisplay -nosplash -nodesktop -r "Parallel_function_AbsDist_Modelv',versionmodel,'_v',versionsubmodel,'_',num2str(i),'"']);
    
    fclose(fileID);
    
end

clear all

filesnotpresent = [];

for filenumber = 1:5000
    
    filename = ['/scratch/ec47/19_10_11_9Mut_SENSITIVITYNEWPRIORS_ABCMv9v1_20000part_64sim_11par_Eps4_Eps2p33_',num2str(filenumber),'.mat'];
    
    if exist(filename, 'file') == 2
        
    else
        
        filesnotpresent = [filesnotpresent,filenumber];
        
    end
    
end

jobstorerun = unique(ceil(filesnotpresent/10))


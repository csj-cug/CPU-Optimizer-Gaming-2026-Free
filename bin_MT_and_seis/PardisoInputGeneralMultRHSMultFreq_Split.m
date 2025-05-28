function [freq_per_comm, nfreq_per_comm] = PardisoInputGeneralMultRHSMultFreq_Split(A,b_mult,complexFlag,symFlag,fileformat,compress,freq,Num_Comm)
%% Write input file for Pardiso
% pardiso expect CSR format for the sparse matrix
% symmetric matrix only stores the upper triangle
% multiple right hand sides are treated as flat vectors b
% i.e. if there are NRHS right hand sides, the vector b size is NRHS*length(b)
% A: input matrix in assembled form (For sparse matrix).
% b: right hand side
% filename: a string for the input file
% complexFlag: 0: A is real, 1: A is complex
% A is symmetry matrix;
% only store the lower triangle part
% symFlag: 1 for symmetry and 0 for general matrix
% freq: a vector stores the frequency value
% ns_per_freq: a vector with the same size of freq, each element indicates
% the number of source for this frequency.

if (size(freq,2))>1
    freq = freq.';
    warning('Frequency should be a column vector instead of a row vector and it is automatically reshaped now.')
end


nfreq = length(freq);


if size(b_mult,1) ~= nfreq
    error('number of cell rows for b_mult should be equal to nfreq');
end

b_mult_transp = b_mult';
clear b_mult;
b_mult = b_mult_transp(:);
clear b_mult_transp;






NRHS=length(b_mult);

if mod(NRHS,nfreq) ~=0
    error('Number of source for each frequency should be an integer');
else
    ns_per_freq = NRHS/nfreq;%number of source for each frequency
end



nfreq_per_comm = zeros(Num_Comm,1);%Number of frequency for each communicator
nfreq_per_comm_round = round(nfreq/Num_Comm);

for icomm = 1:Num_Comm
    nfreq_per_comm(icomm) = nfreq_per_comm_round;
    if icomm == Num_Comm
        nfreq_per_comm(icomm) = nfreq -sum(nfreq_per_comm(1:icomm-1)) ;
    end
end

NRHS_per_comm = nfreq_per_comm*ns_per_freq;


freq_per_comm = cell(Num_Comm,1);%frequency value for each communicator
for icomm = 1:Num_Comm
    if icomm == 1
        freq_per_comm{icomm} = freq(1:nfreq_per_comm(icomm));
    else
        ind = sum(nfreq_per_comm(1:icomm-1));
        freq_per_comm{icomm} = freq(ind+1:ind+nfreq_per_comm(icomm));
    end
    
end


% if sum(ns_per_freq)~=NRHS
%     error('Number of RHS does not match!');
% end


% %write the number of frequency to file
% fid = fopen('Nfreq','w');
% fprintf(fid,'%d\n',nfreq);
% fclose(fid);

%write the frequency
% fid = fopen('FreqAll','w');
% for ifreq = 1:nfreq
%     fprintf(fid,'%10.10g\n',freq(ifreq));
% end
% fclose(fid);

%write number of source for each frequency

% fid = fopen('ns_per_freq','w');
% for ifreq = 1:nfreq
%     fprintf(fid,'%d\n',ns_per_freq(ifreq));
% end
% fclose(fid);



%write number of MPI communicator after split.
fid = fopen('Num_Comm','w');
fprintf(fid,'%d \n',Num_Comm);
fclose(fid);

fid = fopen('FreqFirst','w');
    fprintf(fid,'%10.10g %10.10g\n',[real(freq(1));imag(freq(1))]);
fclose(fid);


b = cell(Num_Comm,1);

for icomm = 1:Num_Comm
    
    if icomm == 1
        for isrc = 1:NRHS_per_comm(icomm)
            b{icomm} = [b{icomm};b_mult{isrc}(:)];
        end
    else
        for isrc = 1:NRHS_per_comm(icomm)
            b{icomm} = [b{icomm};b_mult{isrc+sum(NRHS_per_comm(1:icomm-1))}(:)];
        end
    end
    
    
end



% b=[];
% for isrc = 1:NRHS
%     b = [b;b_mult{isrc}(:)];
% end

[row,col,val]=find(A);

switch symFlag
    case 0
        disp('General matrix is written to input');
        [row,col,val]=sparse_to_csr(A);
        N = length(A);
    case 1
        disp('Symmetric matrix is written to input');
        %        ind = find(row<col);
        %        row(ind)=[];
        %        col(ind)=[];
        %        val(ind)=[];
        AA = triu(A);
        [row,col,val] = sparse_to_csr(AA);
        N = length(AA);
    otherwise
        warning('No such option and general matrix is assumed');
end


NC=length(col);
NR=length(row);

% tic  - WRITE N, NC, NR, NB, NRHS, symm, row, col, AA, b to file
% consider to write binary and/or compress for faster I/O
disp(['Matrix size ' num2str(N) ' with ' num2str(NC) ' nonzero elements']);


% disp('Write FEM information for all frequency together for MPI_COMM_WORLD');
% fid = fopen([filename '_All'],'w');
% fprintf(fid,'%d %d %d %d %d %d %d %d  %d\n',N, NR, NC, NRHS, symFlag, complexFlag, fileformat, compress, nfreq);
% fprintf(fid,'%10.10g %10.10g\n',[real(freq)';imag(freq)']);
% fclose(fid);

disp('Write Input information for all MPI_COMM');
fid = fopen('InputP_All','w');
fprintf(fid,'%d %d %d %d %d %d %d %d  %d\n',N, NR, NC, NRHS, symFlag, complexFlag, fileformat, compress, nfreq);
fprintf(fid,'%10.10g %10.10g\n',[real(freq)';imag(freq)']);
fclose(fid);


disp('Write sparse matrix in CSR format');

if fileformat == 0
    disp('Input matrix is written into text format');
    fid = fopen('InputP_Matrix','w');
    fprintf(fid,'%d\n',row');
    fprintf(fid,'%d\n',col');
    
    if complexFlag == 0
        fprintf(fid,'%10.10g\n',val');
    else
        fprintf(fid,'%10.10g %10.10g\n',[real(val)';imag(val)']);
    end
    fclose(fid);
    
elseif fileformat == 1
    disp('Input matrix is written into HDF5 format');
    fid = H5F.create('InputP_Matrix.h5', 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    % first space - rows
    h5write_int(fid,'ia',row,compress);
    % second space - columns
    h5write_int(fid,'ja',col,compress);
    % third space - matrix
    h5write_cplx(fid,'a',val,compress);
    H5F.close(fid);
else
    error('file format can only be either 0 or 1');
end





for icomm = 1:Num_Comm
    
    filename_comm = ['InputP_Head_Comm_' num2str(icomm-1)];
    fid = fopen(filename_comm,'w');
    fprintf(fid,'%d %d %d %d %d %d %d %d  %d\n',N, NR, NC, NRHS_per_comm(icomm), symFlag, complexFlag, fileformat, compress, nfreq_per_comm(icomm));
    fprintf(fid,'%10.10g %10.10g\n',[real(freq_per_comm{icomm})';imag(freq_per_comm{icomm})']);
    fclose(fid);
    
    if fileformat  == 0
        
        filename_comm = ['InputP_Comm_' num2str(icomm-1)];
        fid = fopen(filename_comm,'w');
        if complexFlag == 0
            fprintf(fid,'%10.10g\n',b{icomm});
        else
            fprintf(fid,'%10.10g %10.10g\n',[real(b{icomm})' ; imag(b{icomm})']);
        end
        fclose(fid);
        
    elseif fileformat == 1
        filename_comm = ['InputP_Comm_' num2str(icomm-1) '.h5'];
        fid = H5F.create(filename_comm, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
        % first space
        h5write_cplx(fid,'b_comm',b{icomm},compress);
        H5F.close(fid);
    else
        error('file format can only be either 0 or 1');
    end
end

1;






% for icomm = 1:Num_Comm
%     filename_comm = ['InputP_Comm_' num2str(icomm-1)];
%     fid = fopen(filename_comm,'w');
%     fprintf(fid,'%d %d %d %d %d %d %d %d  %d\n',N, NR, NC, NRHS_per_comm(icomm), symFlag, complexFlag, fileformat, compress, nfreq_per_comm(icomm));
%     fprintf(fid,'%10.10g %10.10g\n',[real(freq_per_comm{icomm})';imag(freq_per_comm{icomm})']);
%     if complexFlag == 0
%         fprintf(fid,'%10.10g\n',b{icomm});
%     else
%         fprintf(fid,'%10.10g %10.10g\n',[real(b{icomm})' ; imag(b{icomm})']);
%     end
%     
%     fclose(fid);
% end





% if (fileformat == 0)
%     disp('Solver I/O format is text');
%     fid = fopen(filename,'w');
%     fprintf(fid,'%d %d %d %d %d %d %d %d  %d\n',N, NR, NC, NRHS, symFlag, complexFlag, fileformat, compress, nfreq);
%     fprintf(fid,'%10.10g %10.10g\n',[real(freq)';imag(freq)']);
%     fprintf(fid,'%d\n',row');
%     fprintf(fid,'%d\n',col');
%     
%     if complexFlag == 0
%         fprintf(fid,'%10.10g\n',val');
%         fprintf(fid,'%10.10g\n',b);
%     else
%         fprintf(fid,'%10.10g %10.10g\n',[real(val)';imag(val)']);
%         fprintf(fid,'%10.10g %10.10g\n',[real(b)' ; imag(b)']);
%     end
%     fclose(fid);
% elseif (fileformat == 1)
%     disp(['Solver I/O format is HDF5 with compression ' num2str(compress)]);
%     fid = fopen(filename,'w');
%     fprintf(fid,'%d %d %d %d %d %d %d %d  %d\n',N, NR, NC, NRHS, symFlag, complexFlag, fileformat, compress, nfreq);
%     fprintf(fid,'%10.10g %10.10g\n',[real(freq)';imag(freq)']); fclose(fid);
%     
%     
%     
%     
%     fid = H5F.create('Input.h5', 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
%     % first space - rows
%     h5write_int(fid,'ia',row,compress);
%     % second space - columns
%     h5write_int(fid,'ja',col,compress);
%     % third space - matrix
%     h5write_cplx(fid,'a',val,compress);
%     % fourth space - RHS
%     h5write_cplx(fid,'b',b,compress);
%     H5F.close(fid);
% end
% system('sync');
% %system('sleep 10');
% % toc

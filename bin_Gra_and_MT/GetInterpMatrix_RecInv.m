function [Q_Ex, Q_Ey, Q_Ez,Q_Hx,Q_Hy,Q_Hz,Nre]=GetInterpMatrix_RecInv(rec,Ve,t_xyz_mm,Ni,CurlNi,ed_length,TE,inv_flag) %Field Interpolation matrix

nrec = size(rec,1);

ned  = length(ed_length);
Ne   = length(Ve);
Q_Ex=sparse(nrec,ned);
Q_Ey=sparse(nrec,ned);
Q_Ez=sparse(nrec,ned);
Q_Hx=sparse(nrec,ned);
Q_Hy=sparse(nrec,ned);
Q_Hz=sparse(nrec,ned);


for irec = 1:nrec
    [Q_Ex(irec,:), Q_Ey(irec,:), Q_Ez(irec,:), Q_Hx(irec,:), Q_Hy(irec,:), Q_Hz(irec,:)]= ...
        GetInterpMatrix_OneReceiver(rec(irec,1),rec(irec,2),rec(irec,3),Ve,t_xyz_mm,Ni,CurlNi,ed_length,TE);
      
end

switch inv_flag
    case 1 % off diagonal impedance
        Q = [Q_Ex; Q_Ey; Q_Hx; Q_Hy];
        npara = 4;
    case 2 % full impedance
        Q = [Q_Ex; Q_Ey; Q_Hx; Q_Hy];
        npara = 4;
    case 3 % tipper
        Q = [Q_Hx; Q_Hy; Q_Hz];
        npara = 3;
    case 4 % off diagonal impedance + tipper
        Q = [Q_Ex; Q_Ey; Q_Hx; Q_Hy; Q_Hz];
        npara = 5;
    case 5 % full impedance + tipper
        Q = [Q_Ex; Q_Ey; Q_Hx; Q_Hy; Q_Hz];
        npara = 5;
    case 6 % phase tensor
        Q = [Q_Ex; Q_Ey; Q_Hx; Q_Hy];
        npara = 4;
end

fid=fopen('nrec','w');
fprintf(fid,'%d \n',npara*nrec);
fclose(fid);

fid1 = fopen('recnum','w');
fid2 = fopen('ed','w');
fid3 = fopen('Q','w');

%output Q
Nre = 0;
for irec = 1:npara*nrec
    index1 = find(Q(irec,:));
    l1 = length(index1);
    Nre = Nre + l1;
    recnum = irec*ones(size(index1));
    fprintf(fid1,'%d \n',recnum);
    fprintf(fid2,'%d \n',index1);    
    Q1(1,:) = full(Q(irec,:));
    Q2 = [Q1(1,index1)];
    fprintf(fid3,'%10.10g\n',Q2);
end

fclose(fid1);
fclose(fid2);
fclose(fid3);


function [Q_Ex, Q_Ey, Q_Ez, Q_Hx, Q_Hy, Q_Hz]=GetInterpMatrix_OneReceiver(xr,yr,zr,ve,t_xyz_mm,Ni,CurlNi,ed_length,TE)

ned  = length(ed_length);
ind = find(t_xyz_mm(:,1) <= xr & t_xyz_mm(:,2) >= xr & t_xyz_mm(:,3)<=yr&...
        t_xyz_mm(:,4) >= yr & t_xyz_mm(:,5) <= zr & t_xyz_mm(:,6) >= zr);
len_index  =length(ind);

Q_Ex=sparse(1,ned);
Q_Ey=sparse(1,ned);
Q_Ez=sparse(1,ned);
Q_Hx=sparse(1,ned);
Q_Hy=sparse(1,ned);
Q_Hz=sparse(1,ned);
for it=1:12
    Q_Ex1(:,it)=squeeze(Ni{it}(1,14,ind));
    Q_Ey1(:,it)=squeeze(Ni{it}(2,14,ind));
    Q_Ez1(:,it)=squeeze(Ni{it}(3,14,ind));
    Q_Hx1(:,it)=squeeze(CurlNi{it}(1,14,ind));
    Q_Hy1(:,it)=squeeze(CurlNi{it}(2,14,ind));
    Q_Hz1(:,it)=squeeze(CurlNi{it}(3,14,ind));
end

if len_index==1
    %     Q_Hz1=squeeze(CurlNie(index,3,:)).*(TE_Sign(index,:).'); %Sign is
    %     alread taken into account in CurlNie
    Q_Ex(TE(ind,:))=Q_Ex1;
    Q_Ey(TE(ind,:))=Q_Ey1;
    Q_Ez(TE(ind,:))=Q_Ez1;
    Q_Hx(TE(ind,:))=Q_Hx1;
    Q_Hy(TE(ind,:))=Q_Hy1;
    Q_Hz(TE(ind,:))=Q_Hz1;
else %the point is shared by many element
    ve_index = repmat(ve(ind),[1,12]);
    ve_index_sum=repmat(sum(ve_index,1),[len_index,1]);
    
    Q_Ex2=(Q_Ex1.*ve_index)./ve_index_sum;
    Q_Ey2=(Q_Ey1.*ve_index)./ve_index_sum;
    Q_Ez2=(Q_Ez1.*ve_index)./ve_index_sum;
    Q_Hx2=(Q_Hx1.*ve_index)./ve_index_sum;
    Q_Hy2=(Q_Hy1.*ve_index)./ve_index_sum;
    Q_Hz2=(Q_Hz1.*ve_index)./ve_index_sum;
    for it=1:len_index
        Q_Ex(TE(ind(it),:))=Q_Ex(TE(ind(it),:))+Q_Ex2(it,:);
        Q_Ey(TE(ind(it),:))=Q_Ey(TE(ind(it),:))+Q_Ey2(it,:);
        Q_Ez(TE(ind(it),:))=Q_Ez(TE(ind(it),:))+Q_Ez2(it,:);
        Q_Hx(TE(ind(it),:))=Q_Hx(TE(ind(it),:))+Q_Hx2(it,:);
        Q_Hy(TE(ind(it),:))=Q_Hy(TE(ind(it),:))+Q_Hy2(it,:);
        Q_Hz(TE(ind(it),:))=Q_Hz(TE(ind(it),:))+Q_Hz2(it,:);%remember to add them togeter, otherwise, the previous value will be overwritten
    end
end

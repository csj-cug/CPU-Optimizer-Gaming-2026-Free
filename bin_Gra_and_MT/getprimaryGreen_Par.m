function [epx,epy,epz]      = getprimaryGreen_Par(freq,par,ed_center,sl,hl,al,viewflag)
                              
%Get the primary vector and scalr potential for each node
%Paralle version in MATLAB;
%Distributed the edge center on different cores.

global sig_air;

%% Modified on 1/8,2020
%% Considering green3d is not suitable for MT forward problem, the function parameters is modified.
sig_air_top = abs(min(ed_center(:,3)));
hl = [sig_air_top,hl];
sl = [sig_air,sl];
al = [1,al];

ned = size(ed_center,1); %number of edges.
interv = round(ned/8);

x   = cell(8,1);
y   = cell(8,1);
z   = cell(8,1);
E   = cell(8,1);

for it=1:8
    if it<=7
        x{it}       = ed_center((it-1)*interv+1:it*interv,1);
        y{it}       = ed_center((it-1)*interv+1:it*interv,2);
        z{it}       = ed_center((it-1)*interv+1:it*interv,3) + sig_air_top;
    else 
        x{it}       = ed_center((it-1)*interv+1:end,1);
        y{it}       = ed_center((it-1)*interv+1:end,2);
        z{it}       = ed_center((it-1)*interv+1:end,3) + sig_air_top;
    end
    
end

% matlabpool open;
parfor it = 1:8
%     it
    [E{it},~]=green3d(freq,hl,sl,al,x{it},y{it},z{it},par);
end
% matlabpool close;

e  = [];
for it = 1:8
    e = [e;E{it}];
end

% apply MT source
if par(1) == 1 && length(par)==2
    switch par(2)
        case 1
            e(:,[2,3]) =0;
        case 2
            e(:,[1,3]) =0;
        otherwise
            error('par(2) must equal to 1 or 2');
    end
end

epx = e(:,1);
epy = e(:,2);
epz = e(:,3);

if viewflag
    xi = linspace(min(ed_center(:,1)),max(ed_center(:,1)),60);
    yi = linspace(min(ed_center(:,2)),max(ed_center(:,2)),60);
    zi = linspace(min(ed_center(:,3)),max(ed_center(:,3)),60);
    zi = -500:50:1000;
    [xi0, yi0] = meshgrid(xi,yi);
    zi0=0*xi0;
    ind = find(ed_center(:,3)==0);
    [xi, yi, zi] = meshgrid(xi,yi,zi);
    Epxi = griddata(ed_center(:,1),ed_center(:,2),ed_center(:,3),epx,xi,yi,zi);
    Epyi = griddata(ed_center(:,1),ed_center(:,2),ed_center(:,3),epy,xi,yi,zi);
    Epzi = griddata(ed_center(:,1),ed_center(:,2),ed_center(:,3),epz,xi,yi,zi);
    
    Epzi0 = griddata(ed_center(ind,1),ed_center(ind,2),epz(ind),xi0,yi0);
    figure;
    slice(xi,yi,zi,log10(abs(Epxi)),[0],[0],[00]);
    shading flat
    set(gca,'zdir','reverse')
    colorbar;
    title('Epx');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    
    figure;
    slice(xi,yi,zi,log10(abs(Epyi)),[0],[0],[00]);
    shading flat
    set(gca,'zdir','reverse')
    colorbar;
    title('Epy');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    
    figure;
    slice(xi,yi,zi,log10(abs(Epzi)),[0],[0],[00]);
    shading flat
    set(gca,'zdir','reverse')
    colorbar;
    title('Epz');
    xlabel('x');
    ylabel('y');
    zlabel('z');
end




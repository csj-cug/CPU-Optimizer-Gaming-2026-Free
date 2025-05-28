function showdata(intout,intout2)
%This function is written to display the data produced by intem3d.
%if index = 0, show the background field;
%if index = 1, show the anomalous field.
%if index = 3, show the total field.
%if index = 4, show the background field and the anomalous field together
%if index = 5, show the normalized field

if size(intout,2) ~= 10
    error('The data should have 10 columns');
end

comp = unique(intout(:,4));%component
compname = ['Ex';'Ey';'Ez';'Hx';'Hy';'Hz'];
src  = unique(intout(:,5));%source
freq = unique(intout(:,6));%frequency

Ncomp = length(comp);
Nsrc  = length(src);
Nfreq = length(freq);


for ifreq = 1:Nfreq
    indf = find(intout(:,6)==freq(ifreq));
    temp1_intout   = intout(indf,:);
    temp3_intout   = intout2(indf,:);
    for isrc = 1:Nsrc
        inds = find(temp1_intout(:,5)==src(isrc));
        temp2_intout = temp1_intout(inds,:);
        temp4_intout = temp3_intout(inds,:);
        x = unique(temp2_intout(:,1));
        y = unique(temp2_intout(:,2));
        nx = length(x);
        ny = length(y);
        
        for icomp = 1:Ncomp
            indcomp = find(temp2_intout(:,4)==comp(icomp));
                
             
                    %% normalized field
                eval([compname(icomp,:) '_bacreal=temp2_intout(indcomp,7);']);
                eval([compname(icomp,:) '_bacimag=temp2_intout(indcomp,8);']);
                temp1 = eval([compname(icomp,:) '_bacreal;']);
                temp2 = eval([compname(icomp,:) '_bacimag;']);
                realptbac = (reshape(temp1,nx,ny)');
                imagptbac = (reshape(temp2,nx,ny)');
                
                eval([compname(icomp,:) '_anoreal=temp2_intout(indcomp,9);']);
                eval([compname(icomp,:) '_anoimag=temp2_intout(indcomp,10);']);
                temp3 = eval([compname(icomp,:) '_anoreal;']);
                temp4 = eval([compname(icomp,:) '_anoimag;']);
                realptano = (reshape(temp3,nx,ny)');
                imagptano = (reshape(temp4,nx,ny)');
                
                bacfld = realptbac+1i*imagptbac;
                anofld = realptano+1i*imagptano;
                totalfld = bacfld+anofld;
                
                eval([compname(icomp,:) '_bacreal=temp4_intout(indcomp,7);']);
                eval([compname(icomp,:) '_bacimag=temp4_intout(indcomp,8);']);
                temp1 = eval([compname(icomp,:) '_bacreal;']);
                temp2 = eval([compname(icomp,:) '_bacimag;']);
                realptbac = (reshape(temp1,nx,ny)');
                imagptbac = (reshape(temp2,nx,ny)');
                
                eval([compname(icomp,:) '_anoreal=temp4_intout(indcomp,9);']);
                eval([compname(icomp,:) '_anoimag=temp4_intout(indcomp,10);']);
                temp3 = eval([compname(icomp,:) '_anoreal;']);
                temp4 = eval([compname(icomp,:) '_anoimag;']);
                realptano = (reshape(temp3,nx,ny)');
                imagptano = (reshape(temp4,nx,ny)');
                
                bacfld2 = realptbac+1i*imagptbac;
                anofld2 = realptano+1i*imagptano;
                totalfld2 = bacfld2+anofld2;
                
                
                
                fig = figure;
                figure(fig);
                subplot(221);
%                 imagesc(x,y,abs(totalfld)./abs(bacfld));
                imagesc(x,y,real(totalfld)./abs(bacfld));
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Observe NormAmpTotal ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;colormap(jet);
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                subplot(223);
%                 imagesc(x,y,unwrap(unwrap(angle(totalfld))-unwrap(angle(bacfld))));
                imagesc(x,y,imag(totalfld)./abs(bacfld));
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Observe NormalPhaTotal ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;colormap(jet);
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                
                
                
                 subplot(222);
%                 imagesc(x,y,abs(totalfld2)./abs(bacfld2));
                imagesc(x,y,real(totalfld2)./abs(bacfld2));
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Predict NormAmpTotal ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;colormap(jet);
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                subplot(224);
%                 imagesc(x,y,unwrap(unwrap(angle(totalfld2))-unwrap(angle(bacfld2))));
                imagesc(x,y,imag(totalfld2)./abs(bacfld2));
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Precict NormalPhaTotal ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;colormap(jet);
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                
                
                 
            
            
        end
        
    end
end
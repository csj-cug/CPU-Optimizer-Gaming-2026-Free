function showdata(intout,index)
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
    for isrc = 1:Nsrc
        inds = find(temp1_intout(:,5)==src(isrc));
        temp2_intout = temp1_intout(inds,:);
        x = unique(temp2_intout(:,1));
        y = unique(temp2_intout(:,2));
        nx = length(x);
        ny = length(y);
        
        for icomp = 1:Ncomp
            indcomp = find(temp2_intout(:,4)==comp(icomp));
            if index == 0
                eval([compname(icomp,:) '_real=temp2_intout(indcomp,7);']);
                eval([compname(icomp,:) '_imag=temp2_intout(indcomp,8);']);
                temp1 = eval([compname(icomp,:) '_real;']);
                temp2 = eval([compname(icomp,:) '_imag;']);
                realpt = (reshape(temp1,nx,ny));
                imagpt = (reshape(temp2,nx,ny));
                fig = figure;
                figure(fig);
                subplot(211);
                imagesc(x,y,realpt');
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Real Bac ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;colormap (jet)
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                subplot(212);
                imagesc(x,y,imagpt');
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Imag Bac ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
            elseif index == 1
                eval([compname(icomp,:) '_real=temp2_intout(indcomp,9);']);
                eval([compname(icomp,:) '_imag=temp2_intout(indcomp,10);']);
                temp1 = eval([compname(icomp,:) '_real;']);
                temp2 = eval([compname(icomp,:) '_imag;']);
                realpt = (reshape(temp1,nx,ny));
                imagpt = (reshape(temp2,nx,ny));
                fig = figure;
                figure(fig);
                subplot(211);
                imagesc(x,y,realpt');
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Real Ano ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                subplot(212);
                imagesc(x,y,imagpt');
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Imag Ano ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
            elseif index == 3
                eval([compname(icomp,:) '_real=temp2_intout(indcomp,9)+temp2_intout(indcomp,7);']);
                eval([compname(icomp,:) '_imag=temp2_intout(indcomp,10)+temp2_intout(indcomp,8);']);
                temp1 = eval([compname(icomp,:) '_real;']);
                temp2 = eval([compname(icomp,:) '_imag;']);
                realpt = (reshape(temp1,nx,ny));
                imagpt = (reshape(temp2,nx,ny));
                fig = figure;
                figure(fig);
                subplot(211);
                imagesc(x,y,realpt');
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Real Total ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                subplot(212);
                imagesc(x,y,imagpt');
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Imag Total ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                
                
                
            elseif index == 4
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
                

                fig = figure;
                figure(fig);
                subplot(221);
                imagesc(x,y,realptbac);
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Real Bac ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                subplot(223);
                imagesc(x,y,imagptbac);
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Imag Bac ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                
                
                
                subplot(222);
                imagesc(x,y,realptano);
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Real Ano ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                subplot(224);
                imagesc(x,y,imagptano);
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['Imag Ano ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                
                
                
                elseif index == 5
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
                

                fig = figure;
                figure(fig);
                subplot(221);
                imagesc(x,y,abs(anofld)./abs(bacfld));
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['NormAmpAno ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;colormap (jet)
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                subplot(223);
                imagesc(x,y,unwrap(unwrap(angle(anofld))-unwrap(angle(bacfld))));
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['NormalPhaAno ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;colormap (jet)
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                
                
                
                 subplot(222);
                imagesc(x,y,abs(totalfld)./abs(bacfld));
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['NormAmpTotal ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;colormap (jet)
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                subplot(224);
                imagesc(x,y,unwrap(unwrap(angle(totalfld))-unwrap(angle(bacfld))));
                xlabel('x'); ylabel('y'); set(gca,'fontweight','normal');
                title(['NormalPhaTotal ' compname(comp(icomp),:) ' ' num2str(freq(ifreq)) 'Hz Src' num2str(src(isrc))]);
                colorbar; h=colorbar;colormap (jet)
                if comp(icomp)<=3
                    set(get(h,'xlabel'),'string','V/m','fontsize',10);
                else
                    set(get(h,'xlabel'),'string','A/m','fontsize',10);
                end
                axis equal; axis([x(1) x(end) y(1) y(end)]);
                
                
                
            else
                error('The value of index can only be 0 or 1');
            end

            
            
            
        end
        
    end
end
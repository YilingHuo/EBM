function xareaavg= areaavg_lat(x,lat,latlim)
    nsize = size(x);ndim=max(size(nsize));
    if min(nsize)<2
        ndim=ndim-1;
    end
    nlatsize=max(size(size(lat)));
    if min(size(lat))<2
        nlatsize=nlatsize-1;
    end
    latweights=lat*0;
    format long
    latweights(lat>latlim)=cos(pi/180*lat(lat>latlim));
    if (nlatsize<2 & ndim>1)
        latweights =repmat(latweights, 1, nsize(1));
    end
    %Area weight
    if ndim>2
        ntim=nsize(end);
        if ndim>3
            nlev=nsize(3);xareaavg=zeros(nlev,ntim);
            for itim =1:ntim
                for ilev =1:nlev
                    tmp=x(:,:,ilev,itim);ii = ~isnan(tmp);
                    xareaavg(ilev,itim)=sum(tmp(ii).*latweights(ii)) / sum(latweights(ii));
                end
            end
        else
            xareaavg=zeros(ntim);
            for itim =1:ntim
                tmp=x(:,:,itim);ii = ~isnan(tmp);
                xareaavg(itim)=sum(tmp(ii).*latweights(ii)) / sum(latweights(ii));
            end
        end
    else
        ii = ~isnan(x);
        latweights=reshape(latweights,nsize);
        xareaavg=sum(x(ii).*latweights(ii)) / sum(latweights(ii));
    end
end

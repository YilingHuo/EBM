function smean=seasonlonavg(x,m1,m2,optionalremovepoles)
   %%%x shape: lon*lat*time in months;smean shape: lat
    if nargin > 3
      removepoles = optionalremovepoles;
    else
      removepoles = false;
    end
    if m2<m1, nmonth=13+m2-m1;else, nmonth=m2-m1+1;end
    nyear=size(x,3)/12;
    x1=squeeze(mean(x,1));
    smean=squeeze(mean(x1,2))*0;
    if m2<m1
        for iyear =0:nyear-1
            smean=smean+(sum(x1(:,iyear*12+1:iyear*12+m2),2)+sum(x1(:,iyear*12+m1:iyear*12+12),2))/nmonth;
        end
    else
        for iyear =0:nyear-1
            smean=smean+mean(x1(:,iyear*12+m1:iyear*12+m2),2);
        end
    end
    smean= smean/nyear;
    if removepoles
        smean= smean(2:end-1);%%%%remove points at 90S and 90N
    end

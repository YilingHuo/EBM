%%%%% controled input %%%%
load("/pscratch/sd/h/huoyilin/InputsforEBM.mat")
lat=lat(2:end-1);lat_grid = lat(2)-lat(1);
lat_number=length(lat);
dir="/pscratch/sd/h/huoyilin/";
file1="e5.sfc.t2m_ssr_ssrd_tisr_tsr_tsrc_ttr_ttrc.19402022.nc";file2="e5.sfc.d2m_mslhf_msnlwrf_msnlwrfcs_msnswrf_msnswrfcs_msshf_mtnlwrf_mtnlwrfcs_mtnswrf_mtnswrfcs_sp_p62_70_72_74_76.162.19402022.nc";
lat0 = ncread(append(dir,file1),"latitude");
m1=1;m2=12;
seasons='JFMAMJJASOND';
if m2-m1>10
    season='Annual';ssn=1;
elseif m2<m1
    season=strcat(seasons(m1:end),seasons(1:m2));ssn=3;
else
    season=seasons(m1:m2);ssn=2;
end
ifdbk=1;%%% which period of feedbacks and forcing to use; 1 for hist; 2 for ssp245
stryr=1950;endyr=2014;filestryr=1940;
nperiod=2;
SW1   =  ncread(append(dir,file2),"mtnswrf")-ncread(append(dir,file2),"msnswrf")-ncread(append(dir,file2),"msnlwrf")-ncread(append(dir,file2),"mslhf")-ncread(append(dir,file2),"msshf");
Tr=ncread(append(dir,file1),"t2m");
etmp=WVPressure(ncread(append(dir,file2),"d2m"));
hs1=etmp./WVPressure(Tr);
ts_gcm_dhist=ts_gcm_dhist(2:end-1);ts_gcm_dssp245=ts_gcm_dssp245(2:end-1);
smoothwndw=20;
SW(:,1)=interp1(lat0,seasonlonavg(SW1(:,:,12*(stryr-filestryr)+1:12*(endyr+1-filestryr)),m1,m2),lat);
hs(:,1)=interp1(lat0,seasonlonavg(hs1(:,:,12*(stryr-filestryr)+1:12*(endyr+1-filestryr)),m1,m2),lat);
clear SW1,hs1
hs(:,2)=hs(:,1);%%no change in relative humidity
if ifdbk<2
    deltaT=ts_gcm_dhist';
    SW(:,2)=SW(:,1)+Force_dhist(2:end-1)';%SW(:,1)=0;
else
    deltaT=ts_gcm_dssp245';
    SW(:,2)=SW(:,1)+Force_dssp245(2:end-1)';%SW(:,1)=0;
end
To0(1)=areaavg_lat(ts_gcm_hist(2:end-1),lat,-91);To0(2)=To0(1);
tmpp=ncread(append(dir,file2),"sp");
PS(:,1)=interp1(lat0,seasonlonavg(tmpp(:,:,12*(stryr-filestryr)+1:12*(endyr+1-filestryr)),m1,m2),lat);PS(:,2)=PS(:,1);%no change in surface pressure
clear tmpp;

fntsz=12;
figure(2);plot(lat,SW(:,2)-SW(:,1),'LineWidth',3)
xlim([-90 90]);xlabel('latitude','fontsize',fntsz);
title("Forcing difference (W m^{-2})",'fontsize',fntsz);
set(gca,'fontname','times')

Bp(:,1)=-Planckfeedback_hist(2:end-1);Bp(:,2)=-Planckfeedback_ssp245(2:end-1);
%Bp(:,1)=areaavg_lat(Bp(:,1),lat,-91);Bp(:,2)=areaavg_lat(Bp(:,2),lat,-91);%%%No Planck feedback curvature
tmp=Allfeedback_hist-Planckfeedback_hist;TotalF(:,1)=tmp(2:end-1);
tmp=Allfeedback_ssp245-Planckfeedback_ssp245;TotalF(:,2)=tmp(2:end-1);
removeF(:,1)=OHUfeedback_hist(2:end-1);removeF(:,2)=WLfeedback_ssp245(2:end-1);
for iperiod =1:nperiod
    figure(5)
    plot(lat,-Bp(:,iperiod),'LineWidth',3);hold on
    figure(6)
    plot(lat,TotalF(:,iperiod),'LineWidth',3);hold on
end
figure(5);xlim([-90 90]);xlabel('latitude','fontsize',fntsz);
legend('his','ssp245','Location', 'Best','fontsize',fntsz);legend boxoff
title("Planck feedback (W m^{-2} K^{-1})",'fontsize',fntsz);
set(gca,'fontname','times');hold off
figure(6);xlim([-90 90]);xlabel('latitude','fontsize',fntsz);
legend('his','ssp245','Location', 'Best','fontsize',fntsz);legend boxoff
title("Total feedback - Planck (W m^{-2} K^{-1})",'fontsize',fntsz);
set(gca,'fontname','times');hold off

%%%%%%%% defining parameters %%%%%%%%%%%%%
format long
rad  = pi/180;
rlat = (rad.*lat)';
rlat_grid = rad*lat_grid;
del_t= 1. * 3600;%seconds;time step = 3 hours is the longest time step can be used for the lat-dependent Do.
Lv=2.5*10^6;  %latent heat
r=6.378*10^6;   % radius of the earth
Cp=1004;   %heat capacity of air at constant pressure
Rd=287;  %gas constant for dry air 
Rv=461;  %gas constant for water vapor
g=9.8;  % gravity
z=interp1(lat0,squeeze(mean(ncread('/global/cfs/cdirs/m1199/huoyilin/cscratch/ERA5/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc',"Z"),1)),lat);
%%%%cpT + Lq+gz is the MSE at the surface
Rd=287;  %gas constant for dry air 
Rv=461;  %gas constant for water vapor
C=1;  %ocean heat capacity
%% create a lookup table for Ta, given ma after first few steps
Ttb = [200:1.e-1:330]';
for iperiod =1:nperiod
    for ilat =1:lat_number
        matb=Cp*Ttb+Lv*hs(ilat,iperiod)*Rd*WVPressure(Ttb)/(Rv*PS(ilat,iperiod))+z(ilat);
        [l(:,ilat,iperiod)] = polyfit(matb,Ttb, 3);
    end
end
clear Ttb matb
%%%%%%%% end of defining parameters %%%%%%%%%%%%%%%%%

%%%%%%%% a few more terms for the differential equation, not much physical meaning
A=zeros(lat_number,lat_number);
for k=1:lat_number-1
  A(k,k)=-1;
  A(k,k+1)=1;
end
B=zeros(lat_number,lat_number);
for k=2:lat_number
  B(k,k)=1;
  B(k,k-1)=-1;
end

lat_ph=lat+lat_grid/2;
lat_nh=lat-lat_grid/2;
rlat_ph = rad.*lat_ph;
rlat_nh = rad.*lat_nh;
coeff1=(r*rlat_grid)^(-2);
figure(4);clear tmp0
Dofile=["K_search_T42.Q2_T2_PHIS_snowfusion_2xco2_season.mat" "K_search_T42.Q2_T2_PHIS_snowfusion_4xco2_season.mat"];
monthweight=[31 28 31 30 31 30 31 31 30 31 30 31];
dir="/global/cfs/cdirs/m1199/jianlu/fukai/0071_cat_0090/";
file1=["e.e112.E1850C5CN.f19_g16.2xco2/e.e112.E1850C5CN.f19_g16.2xco2.0071_cat_0090.clim.nc" "e.e112.E1850C5CN.f19_g16.4xco2/e.e112.E1850C5CN.f19_g16.4xco2.0071_cat_0090.clim.nc"];
file1=[append(dir,file1(1)) append(dir,file1(2))];
lat1 = ncread(file1(1),"lat");lat1=lat1(2:end-1);
for iperiod =1:nperiod
    load(Dofile(iperiod), "K_search"); 
    tmp=interp1(lat1,K_search(ssn,:),lat);%%%seasonal K_search
    Do(:,iperiod)=tmp'.*PS(:,iperiod)./cos(rlat)/g;Do(find(Do<0))=NaN;
    tmp=squeeze(Do(:,iperiod));
    Do(:,iperiod)=fillmissing(tmp,"nearest");
    tmp0(:,iperiod)=Do(:,iperiod)*g./PS(:,iperiod).*cos(rlat);%.*cos(rlat)
    plot(lat,tmp0(:,iperiod),'LineWidth',3);hold on
end
legend('2xCO2','4xCO2','Location', 'Best','fontsize',fntsz*1.5);legend boxoff ;
ylabel("cos\phi D (m^{2}s^{-1})",'fontsize',fntsz);
xlim([-90 90]);xlabel('latitude','fontsize',fntsz);xticks(-90:30:90)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'fontname','times','fontsize',fntsz*1.5);hold off

%%%%%%% starting the energy balance model
%%%%%Do shape and b decides coverge or not (smaller Do smooth window needs smaller b), Bp decides the shape of T curve
To_0=zeros(lat_number,nperiod)+NaN;
bvalue=-235.3;
b=zeros(lat_number,nperiod)+bvalue;
for iperiod =1:nperiod
    To=To0(iperiod)*ones(lat_number,1);  % initial condition for temperature            
    mo=Cp*To+Lv*hs(:,iperiod)*Rd.*WVPressure(To)./(Rv*PS(:,iperiod))+z';
    for n=1:round(86400*365*10/del_t); % time step is del_t seconds, this is 10 years.
        OLR=Bp(:,ifdbk).*To-b(:,iperiod);  %outgoing longwave; a has now been replaced with Bp--plank feedback coef.
        % general form of temperature tendency due to diffusion (horizontal transport)
        Doperiod=1; %%use the same diffusivity for all iperiods
        %Doperiod=iperiod;
        Diff   =Do(:,Doperiod)*coeff1./cos(rlat).*(cos(rlat_ph)'.*(A*mo)-cos(rlat_nh)'.*(B*mo));
        % for south boundary (first grid point)
        Diff(1)=Do(1,Doperiod)*coeff1/cos(rlat(1))*cos(rlat_ph(1))*(mo(2)-mo(1));
        % for north boundary (last grid point)
        Diff(lat_number)=Do(lat_number,Doperiod)*coeff1/cos(rlat(lat_number))*(-cos(rlat_nh(lat_number))*(mo(lat_number)-mo(lat_number-1)));
        tmp=SW(:,iperiod)-OLR+Diff;
        if iperiod>1
           tmp=tmp+TotalF(:,ifdbk).*deltaT;
           tmp=tmp-removeF(:,ifdbk).*deltaT;
        end
        tmp=tmp./PS(:,iperiod)*g;
        mn=mo+del_t*tmp/C;
        mo=mn;    %updating surface MSE %temperature
        To =  l(1,:,iperiod)'.*mo.^3 + l(2,:,iperiod)'.*mo.^2 +l(3,:,iperiod)'.*mo + l(4,:,iperiod)';
        % break the loop if equilibrium is reached
        convergence(n,iperiod) = max(abs(tmp));
        W=(n>1 && convergence(n,iperiod)<1./10^6);
        if (W==1 | isnan(convergence(n,iperiod)))
            break;
        end
    end
    To_0(:,iperiod)=To;
    % break the loop if nan is found in To
    if(sum(isnan(To))>0 && iperiod<nperiod)
        break
    end
end

figure(1)
for iperiod=2:nperiod
    plot(lat,deltaT,'k','LineWidth',3);hold on
    plot(lat,To_0(:,iperiod)-To_0(:,1),'b','LineWidth',3);
end
xlim([-90 90]);xlabel('latitude','fontsize',fntsz);
if ifdbk<2
    legend('his','EBM','Location', 'Best','fontsize',fntsz);
else
    legend('ssp245','EBM','Location', 'Best','fontsize',fntsz);
end
legend boxoff ;
ylabel("\Delta SAT (K)",'fontsize',fntsz);
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'fontname','times');hold off

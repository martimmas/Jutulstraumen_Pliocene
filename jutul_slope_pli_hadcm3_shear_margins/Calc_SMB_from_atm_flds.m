function as = Calc_SMB_from_atm_flds(UserVar,s,temp_clim,prec_clim,s_ref)
%% SMB calculation 
% - Reads in surface temperature and precipitation monthly climatologies
% - Corrects their values to the actual ice sheet surface based on the temperature lapse rate
% - Assumes all precipitation to fall as snow
% - Computes runoff as 
%       - fraction of precip that falls as liquid (Marsiat, 1994)
%       - melt from the semi-analytical PDD model solution of Calov & Greve (2005; J. Glac)
% - Saves as a griddedInterpolant for later use
%
% pdd model translated from src/general/pdd_m.F90 from SICOPOLIS
% runoff calculation translated from src/ant/boundary_m.F90 from SICOPOLIS
if ~isfield(UserVar.SMB,'Gamma_LR'),gamma_lr = 8e-3; end % lapse-rate temp. correction (8.0 degC/km)
if ~isfield(UserVar.SMB,'Gamma_CC'),gamma_cc = 0.07; end % to be used in the precipitation correction

%temp_clim = temp_clim - 273.15;          % convert from K to deg. C
%prec_clim = prec_clim.*(86400*365*1e-3); % convert from kg m-2 s-1 to m/yr
%prec_clim = prec_clim.*12; % convert m/month to m/yr: used for TraCE and Pliocene exps

% Temperature-lapse-rate-based correction for temperature
corr_fact = -gamma_lr .* (s-s_ref);

% Correcting temperature for elevation
temp_clim(:,1:12) = temp_clim(:,1:12) - corr_fact;
prec_clim(:,1:12) = prec_clim(:,1:12) .* exp(gamma_cc.*corr_fact);      

   
accumulation = nanmean(prec_clim,2);
runoff       = calc_runoff(temp_clim,prec_clim, accumulation);

as = accumulation - runoff;

% figure; PlotMeshScalarVariable(CtrlVar,MUA,as); colorbar;

end
%% Calculate runoff from temperature and precipitation
function runoff = calc_runoff(temp_clim, prec_clim, precipitation)

et = calc_pdd_et(temp_clim);

temp_rain = 7.0;
temp_snow = -10.;
inv_delta_temp_rain_snow = 1.0/(temp_rain-temp_snow);

rho    = 910.;
rhow   = 1000.;
beta1  = 3.0 * (0.001/86400.0)*(rhow/rho);
beta2  = 8.0 * (0.001/86400.0)*(rhow/rho);
Pmax   = 0.6;
% mu     = 9.7155 * (1000.0*86400.0)*(rho/rhow);

melt_star = zeros(size(precipitation));
melt      = zeros(size(precipitation));
runoff    = zeros(size(precipitation));
snowfall  = zeros(size(precipitation));

% Determine fraction of solid precipitation (SOLID_PRECIP==1; Marsiat (1994))
for n=1:12
for i=1:length(temp_clim(:,1))
    if (temp_clim(i,n) >= temp_rain)
        frac_solid = 0.0;
    elseif (temp_clim(i,n) <= temp_snow)
        frac_solid = 1.0;
    else
        frac_solid = (temp_rain-temp_clim(i,n))*inv_delta_temp_rain_snow;
    end
end % for i=1:length(temp_clim(:,1))

snowfall = snowfall + prec_clim(:,n)*frac_solid*(1./12.);
end % for n=1:12

% Determine how much precipitation actually falls as runoff (ABLSURFACE==2)
for i=1:length(precipitation)
if ( precipitation(i) <= (Pmax*snowfall(i)) )

  if ( (precipitation(i)+beta1*et(i)) <= (Pmax*snowfall(i)) ) 
     melt_star(i) = precipitation(i)+beta1*et(i);
     melt(i)      = 0.0;
     runoff(i)    = melt(i);
  else
     melt_star(i) = Pmax*snowfall(i);
     melt(i)      = beta2 ...
                      *(et(i)-(melt_star(i)-precipitation(i))/beta1);
     runoff(i)    = melt(i);
  end % if( (rainfall(i)+beta1*et(i)) <= (Pmax*snowfall(i)) ) 

else

  melt_star(i) = Pmax*snowfall(i);
  melt(i)      = beta2*et(i);
  runoff(i)    = melt(i) + precipitation(i)-Pmax*snowfall(i);

end % if( rainfall(i) <= (Pmax*snowfall(i)) )
end % for i=1:length(temp_clim(:,1))

end % function

%% Calculate excess temperaure (ET) from temp climatology (pdd_m.F90)
function et = calc_pdd_et(temp_clim)
time_year     = 1.0;        % period 1 year
time_year_inv = 1.0/time_year;
d_time        = 1.0/12.0;   % time-step 1 month
inv_sqrt2pi   = 1.0/sqrt(2.0*pi);
s_stat        = 5.0; % from /mnt/3tb/data/sico_in/ant/phys_para_ant.dat
inv_s_stat    = 1.0/5.0;
inv_sqrt2     = 1.0/sqrt(2.0);
pdd_sum       = 0.0;

for n=1:12 % month counter


pdd_sum = pdd_sum ...
          + ( s_stat.*inv_sqrt2pi.*exp(-0.5*(temp_clim(:,n).*inv_s_stat).^2) ...
          + 0.5.*temp_clim(:,n) ...
          .*erfc(-temp_clim(:,n)*inv_s_stat*inv_sqrt2) ) ...
          .*d_time;   % positive degree days (in a * deg C)


end
et = pdd_sum*time_year_inv;   % temperature excess   (deg C)

end
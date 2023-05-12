%
%  D_OBSERVATIONS:  Driver script to create a 4D-Var observations file.
%
%  This a user modifiable script that can be used to prepare ROMS 4D-Var
%  observations NetCDF file.  It sets-up all the necessary parameters and
%  variables. USERS can use this as a prototype for their application.
%

% svn $Id: d_observations.m 1156 2023-02-18 01:44:37Z arango $
%=========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set input/output NetCDF files.

 my_root = '/cluster/projects';

 GRDname = strcat(my_root, '/nn9490k/ROHO800/Grid/ROHO800_grid_fix5.nc');
 OBSname = strcat(my_root, '/nn9297k/shmiak/roho800_data/input_data/da_obs.nc');
 
 % Base date for ROMS observation files: "days since 1968-05-23 00:00:00".
% (The WC13 application has a modified Julian day number as reference
%  time: May-23-1968).

mybasedate = datenum(1948,01,01,0,0,0);

%--------------------------------------------------------------------------
%  Set observation file creation parameter in structure array, S.
%--------------------------------------------------------------------------

%  Observations output file name.

S.ncfile = OBSname;

%  Application grid NetCDF file name.

S.grd_file = GRDname;

%  Set application title.

S.title = 'ROHO800: 800 m resolution';

%  Spherical grid switch.
%
%            [0] Cartesian grid
%            [1] Spherical grid
%  1 is for the analytical grid

S.spherical = 0;

%  Set switches to include the 'obs_lon' and 'obs_lat' variables.
%  This are not used inside ROMS but are needed during pre- and
%  post-processing.

S.do_longitude = 0;
S.do_latitude  = 0;

%  Number of state variables. Usually, ROMS uses 7 variables in the
%  observation state vector (zeta, ubar, vbar, u, v, temperature,
%  and salinity). This number need to be increased when considering
%  additional tracer variables.  This is the value of the NetCDF
%  dimension 'state_variable'.

S.Nstate = 7;

%  Number of data surveys. That is, number of unique survey times
%  in the observational dataset. This is the value of the NetCDF
%  dimension 'survey'.

S.Nsurvey = 1;

%  Total number of observations in space and time. This is the value
%  of the NetCDF dimension 'datum'.  If zero, an unlimited dimension
%  is used.

S.Ndatum = 36500;

%  Set attributes for 'obs_type' variable which assigns the model
%  state variable associated with the observation. Usually, ROMS
%  uses 7 variables in the observation state vector:
%
%      obs_type = 1    free-surface
%      obs_type = 2    vertically integrated u-momentum component
%      obs_type = 3    vertically integrated v-momentum component
%      obs_type = 4    u-momentum component
%      obs_type = 5    v-momentum component
%      obs_type = 6    potential temperature
%      obs_type = 7    salinity
%      obs_type = ...  other passive tracers NAT+1:NT
%  
%  NOTE: We are following the CF compliance rules for variable
%        attributes 'flag_values' and 'flag_meanings'.

S.state_flag_values = 1:1:7;

S.state_flag_meanings=['zeta', blanks(1),                               ...
                       'ubar', blanks(1),                               ...
                       'vbar', blanks(1),                               ...
                       'u', blanks(1),                                  ...
                       'v', blanks(1),                                  ...
                       'temperature', blanks(1),                        ...
                       'salinity'];

%  Set attributes for 'obs_provenance' variable which assigns different
%  flags for each instrument or data source.  This information is used
%  in the observation impact and observation sensitivity analysis. The
%  user has a lot of latitute here.
%
%  NOTE: We are following the CF compliance rules for variable
%        attributes 'flag_values' and 'flag_meanings'. All the
%        blank spaces for each meaning needs to be separate with
%        underscores.

S.origin_flag_values = 1:1:2;

S.origin_flag_meanings = ['ROHO160_temperature', blanks(1),             ...
                          'ROHO160_salinity'];

%  The attribute association between 'flag_values' and 'flag_meanings'
%  is difficult to read when a moderate number of flags are use. To
%  aid the decoding, two readable global attributes are added:
%  'state_variables' and 'obs_provenance' which are stored in
%  S.global_variables and S.global_provenance, respectively.
%
%  NOTE: the 'state_variables' attribute include the units.

newline=sprintf('\n');

S.global_variables=[newline,                                            ...
       '1: free-surface (m) ', newline,                                 ...
       '2: vertically integrated u-momentum component (m/s) ', newline, ...
       '3: vertically integrated v-momentum component (m/s) ', newline, ...
       '4: u-momentum component (m/s) ', newline,                       ...
       '5: v-momentum component (m/s) ', newline,                       ...
       '6: potential temperature (Celsius) ', newline,                  ...
       '7: salinity (nondimensional)'];

S.global_provenance=[newline,                                           ...
       ' 1: Modelled temperature from ROHO160 ', newline,               ...
       ' 2: Modelled salinity from ROHO160'];

%  Set the observation data sources global attribute 'obs_sources'
%  which is stored in S.global_sources (OPTIONAL).

% S.global_sources=[newline,                                              ...
%        'http://opendap.aviso.oceanobs.com/thredds/dodsC ', newline,     ...
%        'http://hadobs.metoffice.com/en3'];

%--------------------------------------------------------------------------
%  Create 4D-Var observations NetCDF file.
%--------------------------------------------------------------------------

% [status]=c_observations(S);
c_observations(S);

%  Update 'units' attribute for time variables. Notice that the time
%  of the observations in the NetCDF file is in DAYS.

avalue=['days since ' datestr(mybasedate,31)];

nc_attadd(OBSname,'units',avalue,'survey_time');
nc_attadd(OBSname,'calendar','gregorian','survey_time');

nc_attadd(OBSname,'units',avalue,'obs_time');
nc_attadd(OBSname,'calendar','gregorian','obs_time');

% [status]=nc_attadd(OBSname,'units',avalue,'survey_time');
% [status]=nc_attadd(OBSname,'calendar','gregorian','survey_time');

% [status]=nc_attadd(OBSname,'units',avalue,'obs_time');
% [status]=nc_attadd(OBSname,'calendar','gregorian','obs_time');


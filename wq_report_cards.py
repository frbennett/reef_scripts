import os
import dsed.const as c
from hydrograph import open_dataset
import pandas as pd
import numpy as np
from glob import glob
import webbrowser

FRAC_TO_PC=100.0

UNITS={
    'Flow':('ML/d',c.CUMECS_TO_ML_DAY,None),
    'Sediment - Fine':('kt',c.KG_TO_KTONS,c.KTONS_PER_ML_TO_MG_PER_L),
    'N_Particulate':('t',c.KG_TO_TONS,c.TONS_PER_ML_TO_MG_PER_L),
    'N_DIN':('t',c.KG_TO_TONS,c.TONS_PER_ML_TO_MG_PER_L),
    'N_DON':('t',c.KG_TO_TONS,c.TONS_PER_ML_TO_MG_PER_L),
    'P_Particulate':('t',c.KG_TO_TONS,c.TONS_PER_ML_TO_MG_PER_L),
    'P_FRP':('t',c.KG_TO_TONS,c.TONS_PER_ML_TO_MG_PER_L),
    'P_DOP':('t',c.KG_TO_TONS,c.TONS_PER_ML_TO_MG_PER_L)
}

CONSTITUENT_RENAMES={
    'TSS':'Sediment - Fine',
    'PN':'P_Particulate',
    'DIN':'N_DIN',
    'DON':'N_DON',
    'PP':'P_Particulate',
    'DIP':'P_FRP',
    'DOP':'P_DOP',
    'Total suspended solids':'Sediment - Fine',
    'Flow_ML':'Flow'
}

REPRESENTIVITY_SCORES={
    'Excellent':4,
    'Good':3,
    'Moderate':2,
    'Indicative':1,
    np.nan:0
}

DEFAULT_WATER_YEARS=[f'{y}-{y+1}' for y in range(1986,2018)]

HG_URL='http://www.flowmatters.com.au/hg_202309'
DASHBOARD='gbr-wq-report-card'

def directories_under(base,full=True):
    full_paths = [d for d in glob(os.path.join(base,'*')) if os.path.isdir(d)]
    if full:
        return full_paths
    return [d.split(os.path.sep)[-1] for d in full_paths]

class WQReportCardBuilder(object):
    def __init__(self,results_dir,data_dest,loads_file,concentrations_file,water_years=DEFAULT_WATER_YEARS,csv_files_have_header=True):
        self.results_dir = results_dir
        self.data_dest = data_dest
        self.csv_files_have_header = csv_files_have_header
        self.water_years=water_years
        self.ds = open_dataset(data_dest,mode='rw')

        self.loads_file = loads_file
        self.concentrations_file = concentrations_file
        self.annual_builder = AnnualComparisonBuilder(results_dir,self.ds,loads_file,concentrations_file,water_years=water_years,csv_files_have_header=csv_files_have_header)

    def clear_dataset(self):
        self.ds.clear()

    def build(self):
        self.clear_dataset()
        self.annual_builder.build()

    def view(self):
        self.ds.host()
        webbrowser.open(f'{HG_URL}/#/{DASHBOARD}?data-config.port={self.ds.port}')

class AnnualComparisonBuilder(object):
    def __init__(self,results_dir,dataset,loads_file,concentrations_file,water_years=DEFAULT_WATER_YEARS,csv_files_have_header=True):
        self.results_dir = results_dir
        self.csv_files_have_header = csv_files_have_header
        self.water_years=water_years
        self.ds = dataset
        self.loads_file = loads_file
        self.concentrations_file = concentrations_file

        self.concentrations = None
        self.loads = None
        self.conc_obs_columns = None
        self.runs = directories_under(results_dir)
        self.region_sites = None

    def initialise_obs(self):
        self.concentrations = pd.read_excel(self.concentrations_file)
        self.conc_obs_columns = [c for c in self.concentrations.columns if '(' in c and not c.endswith('Qty')]

        self.loads = pd.read_csv(self.loads_file)

        regions = set(self.loads.Region).union({'all'})
        concentration_sites = dict(set(zip(self.concentrations['SITE NO'].astype('str'),self.concentrations['SITE NAME'])))
        load_sites = dict(set(zip(self.loads['Node'].astype('str'),self.loads['Site'])))

        self.region_sites = {}
        for r in regions:
            if r =='all':
                local_loads = self.loads
            else:
                local_loads = self.loads[self.loads.Region==r]

            load_sites = dict(set(zip(local_loads['Node'].astype('str'),local_loads['Site'])))
            all_sites = set(concentration_sites.keys()).union(load_sites.keys())
            all_sites = {s:concentration_sites.get(s,load_sites.get(s,None)) for s in all_sites}
            self.region_sites[r] = all_sites

    def site_concentrations(self,site):
        conc_site = self.concentrations[self.concentrations['SITE NO']==site][['DateTime']+self.conc_obs_columns].set_index('DateTime').copy()
        if not len(conc_site):
            return None
        for c in self.conc_obs_columns:
            vals = conc_site[c]
            if vals.dtype=='float':
                continue

            vals = vals.map(lambda v: 0.0 if isinstance(v,str) and v.startswith('<') else v)
            vals = vals.map(lambda v: np.nan if isinstance(v,str) and v.startswith('>') else v)
            
            vals = vals.astype('float')
            conc_site[c] = vals    
        return conc_site

    def site_loads(self,site):
        load_site = self.loads[self.loads.Node==site].set_index('Year')
        if not len(load_site):
            return None

    #     load_site['RepresentivityScore'] = load_site['RepresentivityRating'].map(REPRESENTIVITY_SCORES)
    #     load_site = load_site.sort_values('RepresentivityScore',ascending=True)
    #     load_site = load_site.reset_index().drop_duplicates(subset=['Year']).set_index('Year')
        columns = set(load_site.columns) - {'Region','Site','Node','RepresentivityRating','RepresentivityScore'}
        return load_site[columns]

    def load_csv(self,fn):
        if self.csv_files_have_header:
            df = pd.read_csv(fn, index_col = 'Date', parse_dates=True)
        else:
            df = pd.read_csv(fn,index_col=0,header=None,parse_dates=True)
            df = df.rename(columns={1:fn.split(os.path.sep)[-1].replace('.csv','')})
        return df

    def load_obs(self):
        self.ds.rewrite(False)
        all_sites = self.region_sites['all']
        for site,site_label in all_sites.items():
            conc_obs = self.site_concentrations(site)
            if conc_obs is None: continue

            for col in conc_obs.columns:
                ts = conc_obs[col].dropna()
                ts.name = 'data'
                if not len(ts):
                    continue
                constituent = col.split(' (')[0]#.replace('°C','degC').replace('µ','u')
                constituent = CONSTITUENT_RENAMES.get(constituent,constituent)
                units = col.split(' (')[-1].replace(')','')
                self.ds.add_timeseries(ts,site=site,label=site_label,constituent=constituent,units=units,sourced='observation',kind='concentration')

        for site,site_label in all_sites.items():
            load_obs = self.site_loads(site)
            if load_obs is None: continue
            for col in load_obs.columns:
                df = load_obs[[col]].copy()
                if col=='TSS':
                    df[col] *= c.TONS_TO_KTONS
                df = df.rename(columns={col:'observed-load'})
                constituent = CONSTITUENT_RENAMES.get(col,col)
                df = df.reindex(self.water_years)
                self.ds.add_table(df,site=site,label=site_label,constituent=constituent,sourced='observation',kind='load')
        self.ds.rewrite(True)

    def load_run(self,r,sites,**tags):
        ts_dir = os.path.join(r,'TimeSeries')
        if not os.path.exists(ts_dir):
            print(f'No time series for {r}')
            return

        constituents = directories_under(ts_dir,full=False)

        sites_seen = set()

        for site_name, site_label in sites.items():
            for constituent in constituents:
                con_dir = os.path.join(ts_dir,constituent)
                if constituent=='Flows':
                    constituent='Flow'
                ts_pattern = f'{constituent}_{site_name}_*.csv'
                files = list(glob(os.path.join(con_dir,ts_pattern)))
                if len(files)==0:
                    ts_pattern = f'{constituent}_GS{site_name}_*.csv'
                    files = list(glob(os.path.join(con_dir,ts_pattern)))
                    if len(files)==0:
                        continue

                sites_seen = sites_seen.union({(site_name,site_label)})
                if len(files)!=1:
                    print(site_name,constituent,files)
                    assert False

                fn = files[0]
                try:
                    df = self.load_csv(fn)
                except:
                    print(f'Error reading {fn}')
                    raise
                ts = df[df.columns[0]]
                if constituent in UNITS:
                    new_units,load_conversion,_ = UNITS[constituent]
                    try:
                        ts = ts.copy() * load_conversion
                    except:
                        print(f'Error converting {ts.name} to {new_units}')
                        raise
                    ts.name = f'{constituent}_{site_name}_{new_units}'
                    
                    units = new_units
                else:
                    units = ts.name.split('_')[-1]

                ts_tags = dict(site=site_name,label=site_label,constituent=constituent,units=units,sourced='model',kind='load',**tags)
                self.ds.add_timeseries(ts,timestep='day',**ts_tags)
                monthly = ts.resample('1M').sum()
                self.ds.add_timeseries(ts,timestep='month',**ts_tags)
                
                wy = ts.index.map(lambda dt: f'{dt.year}-{dt.year+1}' if dt.month > 6 else f'{dt.year-1}-{dt.year}')
                wy_ts = ts.groupby(wy).sum()
                wy_ts.name = 'modelled-load'
                tbl_len_b4 = len(self.ds.index['tables'])
                self.ds.add_table(wy_ts,timestep='wateryear',**ts_tags)
                assert len(self.ds.index['tables']) == (tbl_len_b4+1)
        return sites_seen

    def load_runs(self):
        self.ds.rewrite(False)
        for run in self.runs:
            run_name = run.split(os.path.sep)[-1]
            region = run_name.split('_')[0]
            tags = {
                'run':run_name,
                'region':region
            }
            print(f'{region}: Loading run {run}')
            sites_seen = self.load_run(run,self.region_sites[region],**tags)
            sites_seen = pd.DataFrame(sites_seen,columns=['Site','Name'])
            self.ds.add_table(sites_seen,purpose='site-lookup',**tags)
        
        self.ds.rewrite(True)

    def summarise_site(self,s):
        def make_df(columns):
            columns = [(con,series) for con,series in columns if series is not None and len(series)]
            if not len(columns):
                return None
            return pd.DataFrame(dict(columns))

        def add_table(columns,**tags):
            df = make_df(columns)
            if df is None:
                return
            self.ds.add_table(df,**tags)

        table_constituents = ['Flow'] + list(self.ds.tag_values('constituent',site=s)-{'Flow','summary','error'})
        flow_obs=None
        flow_mod=None
        mod_load_columns=[]
        err_load_columns=[]

        mod_conc_columns=[]
        err_conc_columns=[]

        for con in table_constituents:
            mod = self.ds.get_tables(constituent=con,site=s,sourced='model')
            if len(mod)==1:
                mod = mod[0]['modelled-load']
            else:
    #             print(f'No modelled {con} in {s}')
                continue
            mod_load_columns.append((con,mod))

            obs = self.ds.get_tables(constituent=con,site=s,sourced='observation',kind='load')
            if len(obs)==1:
                obs = obs[0]['observed-load']
                error = FRAC_TO_PC*(mod-obs)/obs
                err_load_columns.append((con,error))
            else:
                obs = None
                err_load_columns.append((con,None))

            if con=='Flow':
                flow_obs = obs
                flow_mod = mod
            else:
                concentration_conversion = c.KG_PER_ML_TO_MG_PER_L
                if con in UNITS:
                    concentration_conversion = UNITS[con][2]

                mod_conc = mod / flow_mod
                mod_conc *= concentration_conversion
                mod_conc_columns.append((con,mod_conc))

                if obs is None:
                    err_conc_columns.append((con,None))
                else:
                    obs_conc = obs / flow_obs
                    obs_conc *= concentration_conversion
                    conc_error = FRAC_TO_PC*(mod_conc-obs_conc)/obs_conc
                    err_conc_columns.append((con,conc_error))

        add_table(mod_load_columns,site=s,sourced='model',constituent='summary',timestep='wateryear',kind='load')
        add_table(err_load_columns,site=s,sourced='model',constituent='error',timestep='wateryear',kind='load')
        add_table(mod_conc_columns,site=s,sourced='model',constituent='summary',timestep='wateryear',kind='concentration')
        add_table(err_conc_columns,site=s,sourced='model',constituent='error',timestep='wateryear',kind='concentration')

    def summarise_sites(self):
        sites=list(self.ds.tag_values('site'))
        for site in sites:
        #     print(site)
            self.summarise_site(site)


    def build(self):
        self.initialise_obs()
        self.load_obs()
        self.load_runs()
        self.summarise_sites()

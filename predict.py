# Basic imports
import numpy as np
import datetime as dt
import astropy.units as u
from warnings import warn
# Astropy imports
from astropy.time import Time
from astropy.utils.exceptions import AstropyUserWarning, AstropyWarning
from astropy.table import vstack, Table, Column
# Import astroplan tools
from astroplan import (FixedTarget, Observer, EclipsingSystem,
                       is_event_observable, PeriodicEvent)

# Import astroplan constraints
# Other constraints are available
from astroplan import (PrimaryEclipseConstraint, 
                       AtNightConstraint, AltitudeConstraint,
                       TimeConstraint, AirmassConstraint, PhaseConstraint)

# Local imports
from list_of_constraints import List_of_constraints 
from tools import min_start_times, max_end_times


# Function

def load_from_masterfile(*args):
    
    # Need to be installed, so import if necessary
    from masterfile.archive import MasterFile
    
    # Load
    return MasterFile.load().by_pl_name(*args)


# Classes

class Prediction():

    def __init__(self, time_range, targets, site='cfht', constraints=None, supp_cols=None):
        
        
        # Get infos from the MasterFile
        if isinstance(targets, list):
            info = load_from_masterfile(*targets)
            supp_cols = supp_cols or ['pl_orbper','st_j','st_h',
                                      'ra','dec','pl_eqt','st_teff']
        else:
            info = targets.copy()
            if supp_cols is None:
                supp_cols = list(info.keys())
                supp_cols.remove('pl_name')
            
        # Default constraint
        if constraints is None:
            constraints = [AtNightConstraint.twilight_nautical(),
                           AirmassConstraint(max=2.5)]

        # Define time constraint and append it
        t1, t2 = Time(time_range)
        constraints.append(TimeConstraint(t1,t2))
        
        # Convert to List_of_constraints (useful to print and save)
        constraints_list = List_of_constraints(constraints)
        
        # Save infos
        self.info = info
        self.constraints = constraints_list
        self.meta = {'Time_limits': [t1, t2],
                     'Target_list': info['pl_name'].tolist(),
                     'Site': site,
                     **constraints_list.show()}
        self.supp_cols = supp_cols
        self.info_cols = ['pl_name'] + supp_cols
        self.obs = Observer.at_site(site)
#         self.n_eclipses = n_eclipses

        # Resolve targets
        self.targets = [self.resolve_target(i)
                        for i in range(len(targets))]
        
    def resolve_target(self, itar):
        
        info = self.info[itar]
        
        # First simply try to resolve with the name and return
        try:
            target = FixedTarget.from_name(info['pl_name'])
            return target
        except NameResolveError as e:
            print(e)
            pass
        
        # If failed, try to change the name
        try:
            try_name = ' '.join(info['pl_name'].split(' ')[:-2])
            print("Trying with {}".format(try_name))
            target = FixedTarget.from_name(try_name)
        # Finally, try with ra and dec
        except NameResolveError as e:
            print(e)
            print('Searching with RA and dec')
            ra = info['ra'].quantity
            dec = info['dec'].quantity
            coord = SkyCoord(ra=ra, dec=dec)
            target = FixedTarget(coord=coord,
                                 name=info['pl_name'])
            
        return target
    
    def add_column(self, value, key, description, unit=None, message=None):

        info = self.info
        
        if message is None:
            message = "'{}' not available for all systems."  \
                    + " Please specify a value as input."
            message.format(key)
        
        if value is not None:
            try:
                info[key] = Column(value,
                                   unit=value.unit,
                                   description=description)
            except AttributeError:
                info[key] = Column(value,
                                   unit=unit,
                                   description=description)
        elif info[key].mask.any():
            raise ValueError(message)

    def __repr__(self):
        '''
        Represent as a table
        '''
        out = 'Event Prediction \n'
        out += '-------------------\n'
        out += '\n'.join([key + ': ' + str(self.meta[key])
                          for key in self.meta.keys()]) + '\n'
        out += '-------------------\n'
        table = Table(self.info[self.info_cols], meta=self.meta)
        out += table.__repr__()

        return out
        
        
class PredictPhase(Prediction):
    '''
        args :
        -----
        time_range: list of dates
        targets: astropy.table like object of target information.
                 Can also take a list of string of the names of the targets.
        kwargs : 
        --------
        phase_zero: list or quantity object
            list of the time (bjd) of phases zero.
            If not specified, take transit mid point
        site: str, default is 'cfht'
            observing site
        constraints: list, default is None
            observing constraints
        supp_cols: list, default is None  
            supplementary columns to include in the output table.
        '''
    
    def __init__(self, *args, phase_zero=None, **kwargs):
        
        super().__init__(*args, **kwargs)
        info = self.info
        
        # Add phase zero column
        key = 'pl_phase_zero'
        description = 'equivalent to mid-transit time for'  \
                    + 'a non eclipsing system. Use with caution'
        
        if phase_zero is not None:
            try:
                info[key] = Column(phase_zero,
                                   unit=phase_zero.unit,
                                   description = description)
            except AttributeError:
                info[key] = Column(phase_zero,
                                   unit=u.d,
                                   description = description)
        else:
            try:
                info[key]
            except KeyError:
                if info['pl_tranmid'].mask.any():
                    raise ValueError(
                        "'pl_tranmid' not available for all systems." +
                        " Please specify a value for input 'phase_zero'.")
                warn("'{}' not available.".format(key)
                     + " Taking 'pl_tranmid' instead.")
                info[key] = Column(info['pl_tranmid'].quantity,
                                   description = description)
                
        # Save other attributes
        info_cols = ['pl_name', 'pl_phase_zero'] + self.supp_cols
        self.info_cols = list(set(info_cols))  # Make sure cols are unique
        

    def predict(self, phase_range, obs_time=3.*u.h, dt_grid=0.5*u.h):
        '''
        Parameters:
        -----------
        phase_range: list, len=2
            Phase range
        obs_time: quantity or float
            minimum required observing time
        dt_grid: quantity or float
            grid steps used to compute observability between phase range.
        '''
        if not hasattr(obs_time, 'unit'):
            obs_time = obs_time * u.h
        if not hasattr(dt_grid, 'unit'):
            dt_grid = dt_grid * u.h

        t1, t2 = self.meta['Time_limits']
        info = self.info
        n_eclipses = 500  # TODO: COuld be computed according to period and time range
        constraints_list = self.constraints
        obs = self.obs
        supp_cols = self.supp_cols
            
        # Define needed quantities based on planets infos
        # Must be quatities arrays (astropy)
        # Here we use a given astropy Table (info) to get the infos

        epoch, period = [info[k_col].quantity for k_col 
                         in ('pl_phase_zero', 'pl_orbper')]
        epoch = Time(epoch, format='jd')
        pl_name = info['pl_name']
        
        window_start = t1
        
        # Init output table
        col_names = ('pl_name',
                     'Obs_start',
                     'Phase_start',
                     'Obs_end',
                     'Phase_end',
                     'mid_tr',
                     'AM_mid_tr',
                     'moon',
                     *supp_cols
                    )
        
        meta = {'Phase_range': phase_range,
                **self.meta
               }
        full_table = Table()
        
        # Iterations on the targets
        for itar, target in enumerate(self.targets):

            # -------------------------
            # Steps to predict transits
            # -------------------------

            d_phase = (dt_grid /period[itar]).decompose().value
            [p1, p2] = phase_range
            p1 += d_phase/10  # Make sure it's not on the boundary
            p2 -= d_phase/10
            phase_grid = np.arange(p1, p2, d_phase)
            phase_grid = phase_grid*period[itar] + epoch[itar]

            while True:
                # Find all events for each phase point in grid ...
                t_grid = []
                for phase in phase_grid:
                    # Define a system for each phase point
                    sys = EclipsingSystem(primary_eclipse_time=phase,
                                          orbital_period=period[itar])
                    # Compute all events and save
                    t_temp = sys.next_primary_eclipse_time(window_start,
                                                           n_eclipses=n_eclipses)
                    t_grid.append(t_temp.jd)
                # Convert to Time object
                t_grid = Time(t_grid, format='jd').T
                
                # Do so until t2 is passed
                if t_grid[-1,-1] > t2: break
                # or add eclipses to pass t2 and recompute t_grid
                else: n_eclipses += 500

            if (np.diff(t_grid.jd, axis=-1) <= 0).any():
                message = 'Time limit t1 falls into phase range.'  \
                        + ' This will be corrected eventually.'  \
                        + ' For now, please change the time limit t1.'
                raise ValueError(message)

            t_grid = t_grid[(t_grid < t2).any(axis=1)]

            events = []
            for grid in t_grid:
                index = is_event_observable(constraints_list, obs, target, times=grid).squeeze()
                if index.any():
                    events.append(np.mean(grid[index].jd))
            events = Time(events, format='jd')

            # Finally add phase constraint
            sys = PeriodicEvent(epoch=epoch[itar], period=period[itar])

            final_constraints = [*constraints_list,
                PhaseConstraint(sys, *phase_range)]

            # TODO: Add something to check the dt in min_start_times (can bug if dt_grid too small)
            obs_start = min_start_times(final_constraints, obs, target, events)
            obs_end = max_end_times(final_constraints, obs, target, events)
            baseline = obs_end - obs_start
            t_mid = obs_start + baseline/2

            index = (obs_end - obs_start) > obs_time

            # -------------------
            # End of steps to predict events
            # -------------------

            # Put the infos in a table and stack it to the full table
            if index.any():
                name = np.repeat(target.name,index.sum()).astype(str)
                moon = obs.moon_illumination(t_mid[index])
                phase_start = sys.phase(obs_start[index])
                phase_end = sys.phase(obs_end[index])
                AM_mid = obs.altaz(t_mid[index], target).secz
                supp = [np.repeat(info[key][itar],index.sum())
                        for key in supp_cols]
                cols = [name,
                        obs_start[index].iso,
                        phase_start,
                        obs_end[index].iso,
                        phase_end,
                        t_mid[index].iso,
                        AM_mid,
                        moon,
                        *supp
                       ]
                table_sys = Table(cols, names=col_names, masked=True)
                full_table = vstack([table_sys, full_table])
            else:
                warn('No event found for ' + sys.name, AstropyUserWarning)

        if full_table:
            full_table.sort('mid_tr')
            full_table.meta = meta
        else:
            warn('No event found at all', AstropyUserWarning)
            
        return full_table
    
class PredictTransit(Prediction):
    
    def __init__(self, *args, t0=None, duration=None, **kwargs):
    
        super().__init__(*args, **kwargs)
        info = self.info
        
        # Add mid transit time column
        key, description = 'pl_tranmid', 'mid-transit time'
        self.add_column(t0, key, description, unit=u.d)
        
        # Add transit duration column
        key, description = 'pl_trandur', 'transit duration'
        self.add_column(duration, key, description, unit=u.h)
        
#         if t0 is not None:
#             try:
#                 info[key] = Column(t0,
#                                    unit=t0.unit,
#                                    description = description)
#             except AttributeError:
#                 info[key] = Column(t0,
#                                    unit=u.d,
#                                    description = description)
#         elif info[key].mask.any():
#             raise ValueError(
#                 "'pl_tranmid' not available for all systems." +
#                 " Please specify a value for input 't0'.")
                
        # Save other attributes
        info_cols = ['pl_name', 'pl_tranmid'] + self.supp_cols
        self.info_cols = list(set(info_cols))  # Make sure cols are unique
        

    def predict(self):
        
#         # Make sure baseline has units
#         if baseline is not None:
#             try:
#                 baseline.unit
#             except AttributeError:
#                 warn("No units specified for input 'baseline'."
#                      +" Assuming hours.")
#                 baseline = baseline * u.h
            
        
        # Inputs from object's attributes
        t1, t2 = self.meta['Time_limits']
        info = self.info
        n_eclipses = self.n_eclipses
        constraints_list = self.constraints
        obs = self.obs
        supp_cols = self.supp_cols
            
        # Define needed quantities based on planets infos
        # Must be quatities arrays (astropy)
        # Here we use a given astropy Table (info) to get the infos

        epoch, period, transit_duration =   \
            [info[k_col].quantity for k_col 
             in ('pl_tranmid', 'pl_orbper', 'pl_trandur')]
        epoch = Time(epoch, format='jd')
        pl_name = info['pl_name']
        
        observing_time = t1
        
        
        # Init output table
        col_names = ('pl_name',
                     'mid_tr',
                     'AM_mid_tr',
                     'tr_start',
                     'tr_end',
                     'AM_tr_start',
                     'AM_tr_end',
                     'start',
                     'end',
                     'AM_start',
                     'AM_end',
                     'Obs_start',
                     'Baseline_before',
                     'Obs_end',
                     'Baseline_after',
                     'moon',
                     *supp_cols
                    )
        meta = {
                **self.meta
               }
        full_table = Table()
        
        # Iterations on the targets
        for itar, target in enumerate(self.targets):

            # -------------------------
            # Steps to predict transits
            # -------------------------
            
            # Define system
            sys = EclipsingSystem(primary_eclipse_time=epoch[itar],
                                  orbital_period=period[itar],
                                  duration=transit_duration[itar],
                                  name=target.name
                                 )
            
            # Find all events ...
            while True:
                t_mid = sys.next_primary_eclipse_time(observing_time, n_eclipses=n_eclipses)
                # ... until t2 is passed
                if t_mid[-1] > t2: break
                # or add eclipses to pass t2
                else: n_eclipse += 500

            # Remove events after time window
            t_mid = t_mid[t_mid < t2]
            
            # Number of events
            n_event, = t_mid.shape
            
            # Get ingress and egress times
            t1_t4 = sys.next_primary_ingress_egress_time(observing_time, n_eclipses=n_event)

            # Which mid transit times are observable
            i_mid = is_event_observable(constraints_list, obs, target,
                                        times=t_mid).squeeze()
            
            # Which ingress are observable ...
            i_t1 = is_event_observable(constraints_list, obs, target,
                                        times=t1_t4[:,0]).squeeze()
            # ... when mid transit is not.
            i_t1 = i_t1 & ~i_mid
            
            # Which egress are observable ...
            i_t4 = is_event_observable(constraints_list, obs, target,
                                        times=t1_t4[:,1]).squeeze()
            # ... when mid transit and ingress is not.
            i_t4 = i_t4 & ~i_mid & ~i_t1
            
            # Keep events where ingress, mid_transit or egress is observable
            index = i_mid | i_t1 | i_t4
            
            # Get observability for these events.
            # Starting point to compute the observability
            t_obs = np.concatenate([t_mid[i_mid], t1_t4[i_t1,0], t1_t4[i_t4,1]])
            t_obs = Time(t_obs).sort()
            # Get observability range for each of these events
            obs_start = min_start_times(constraints_list, obs, target, t_obs)
            obs_end = max_end_times(constraints_list, obs, target, t_obs)

            # -------------------
            # End of steps to predict transits
            # -------------------

            # Put the infos in a table and stack it to the full table
            if index.any():
                name = np.repeat(sys.name,index.sum()).astype(str)
                moon = obs.moon_illumination(t_mid[index])
                AM_mid = obs.altaz(t_mid[index], target).secz
                AM_t1_t4 = obs.altaz(t1_t4[index], target).secz
        #         AM_base = obs.altaz(t_baseline[index], target).secz
                obs_start = min_start_times(constraints_list, obs, target, t_mid[index])
                baseline_before = (t1_t4[index][:,0] - obs_start).to('min')
                obs_end = max_end_times(constraints_list, obs, target, t_mid[index])
                baseline_after = (obs_end - t1_t4[index][:,1]).to('min')
                supp = [np.repeat(data[key][itar],index.sum())
                        for key in supp_cols]
                cols = [name,
                        t_mid[index].iso,
                        AM_mid,
                        *t1_t4[index].T.iso,
                        *AM_t1_t4.T,
#                         *t_baseline[index].T.iso,
                        *AM_base.T,
                        obs_start.iso,
                        baseline_before,
                        obs_end.iso,
                        baseline_after,
                        moon,
                        *supp
                       ]
                table_sys = Table(cols, names=col_names, masked=True)
                full_table = vstack([table_sys, full_table])
            else:
                warnings.warn('No event found for '+sys.name, AstropyUserWarning)

        if full_table:
            full_table.sort('mid_tr')
            full_table.meta = meta
        else:
            warnings.warn('No event found at all', AstropyUserWarning)

        return full_table

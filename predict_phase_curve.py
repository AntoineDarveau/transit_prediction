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

    def __init__(self, time_range, targets, site='cfht', constraints=None, supp_cols=None, n_eclipses=1000):
        
        
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
        self.n_eclipses = n_eclipses

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
        
        t1, t2 = self.meta['Time_limits']
        info = self.info
        n_eclipses = self.n_eclipses
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
        
        observing_time = t1
        
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

            dt_grid = 0.5*u.h
            d_phase = (dt_grid /period[itar]).decompose().value
            [p1, p2] = phase_range
            p1 += d_phase/10  # Make sure it's not on the boundary
            p2 -= d_phase/10
            phase_grid = np.arange(p1, p2, d_phase)
            phase_grid = phase_grid*period[itar] + epoch[itar]

            t_grid = []
            for phase in phase_grid:
                sys = EclipsingSystem(primary_eclipse_time=phase,
                                      orbital_period=period[itar])
                t_temp = sys.next_primary_eclipse_time(observing_time,
                                                       n_eclipses=n_eclipses)
                t_grid.append(t_temp.jd)
            t_grid = Time(t_grid, format='jd').T


            if t_grid[-1,-1] < t2:
                warn('end time ('+ t2.value +
                     ') is passed the last computed event time (' +
                      t_mid[-1].value+')\n' +
                     '\t You can change the n_eclipse kwarg ' +
                     'value or choose a different window (start or end time)',
                     AstropyUserWarning
                    )

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

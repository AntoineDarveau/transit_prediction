from astropy.time import Time
import astropy.units as u
import numpy as np
from astroplan import is_event_observable


def min_start_time(constraints_list, obs, target, time, dt=1*u.min):
    
    time = Time(time)
    evening = obs.sun_set_time(time, which='previous')

    time_grid = np.arange(evening.jd, time.jd, dt.to('day').value) * u.d
    time_grid = Time(time_grid, format='jd')
    
    index = is_event_observable(constraints_list,
                                obs, target,
                                times=time_grid).squeeze()
    
    return time_grid[index][0]

def max_end_time(constraints_list, obs, target, time, dt=1*u.min):
    
    time = Time(time)
    morning = obs.sun_rise_time(time, which='next')

    time_grid = np.arange(time.jd, morning.jd, dt.to('day').value) * u.d
    time_grid = Time(time_grid, format='jd')
    
    index = is_event_observable(constraints_list,
                                obs, target,
                                times=time_grid).squeeze()
    
    return time_grid[index][-1]

def min_start_times(constraints_list, obs, target, times, dt=1*u.min):
    
    out = []
    for time in times:
        out.append(min_start_time(constraints_list, obs, target, time, dt=1*u.min))
        
    return Time(out)

def max_end_times(constraints_list, obs, target, times, dt=1*u.min):
    
    out = []
    for time in times:
        out.append(max_end_time(constraints_list, obs, target, time, dt=1*u.min))
        
    return Time(out)
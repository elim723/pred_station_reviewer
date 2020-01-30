#!C:/Users/elim.thompson/AppData/Local/Programs/Python/Python37/python

####
#### By Elim Thompson (01/29/2020)
####
#### The python script uses reviewer classes to check for invalid station metadata that are
#### required for prediction calculation. This script can be re-used regularly to monitor
#### changes on station metadata API.
####  
#### For prediction calculation, we need the following information:
####   1. (tides & current harmonic / hydraulic) time zones 
####   2. (tides & current harmonic / hydraulic) harmonic constituents
####   3. (tides & current subordinate) valid reference stations
####   4. (tides & current subordinate) prediction offsets
####   5. (tides harmonic) datums
####   6. (current harmonic / hydraulic) major axis directions
##############################################################################################

import reviewer
from reviewer import tide_reviewer, current_reviewer

## A short lambda function to reformat invalid ref info outputs from a ref-sub map
## to a string 'ref: sub1, sub2' for printing.
reformat_ref_info = lambda ref_info: ['{0}: {1}'.format (ref, ', '.join (slist))
                                      for ref, slist in ref_info.items()]

def print_invalid_stations (stations, info_name, station_type):
    ''' Function to print a list of stations with invalid information. If input list is empty,
        a 'all good' message will be printed. For invalid ref stations, we print the station
        list based on ref in a format of 'ref: sub1, sub2' where sub1 and sub2 are the
        subordinate stations that reference to ref.

        input param
        -----------
        stations   (array): a list of stations with invalid information
        info_name    (str): the name of the invalid information e.g. time zones / datums / etc
        station_type (str): type of stations e.g. harmonic tide / subordinate current / etc
    '''
    ## a dictionary to contain all elements of the print message
    params = {'n_stations'  :len (stations),
              'station_type':station_type  ,
              'info_name'   :info_name      }

    ## In case of invalid ref info,
    ##  1. convert the input 'stations' from dictionary to array of 'ref: sub1, sub2'.
    ##  2. count how many subordinate stations (in values) in total
    if 'ref' in info_name:
        params['n_stations'] = sum ([len (slist) for ref, slist in stations.items()])
        stations = reformat_ref_info (stations)

    if params['n_stations'] == 0:
        print ('    All {station_type} stations have valid {info_name}!'.format (**params))
        return

    ## Print messages on console.
    message  = '    +---------------------------------------------------------------------------------\n'
    message += '    | {n_stations} {station_type} stations with invalid {info_name}:\n'.format (**params)
    for station in stations:
        message += '    |     -- {0}\n'.format (station)
    message += '    +---------------------------------------------------------------------------------'
    print (message)

if __name__ == '__main__':

    ## Create reviewers for tide & current prediction stations
    print ('1. Create reviewers for tides & current prediction stations')
    tider, currenter = tide_reviewer(), current_reviewer()
    print ('')
    
    ## Check time zones for tides & current harmonic & hydraulic stations
    print ('2. Check for time zones - required for harmonic & hydraulic stations')
    invalid_tide = tider.get_stations_with_invalid_timezones()
    print_invalid_stations (invalid_tide, 'time zones', 'harmonic tide')
    print ('')
    invalid_current = currenter.get_stations_with_invalid_timezones()
    print_invalid_stations (invalid_current, 'time zones', 'harmonic & hydraulic current')
    print ('')

    ## Check harmonic constituents for tides & current harmonic & hydraulic stations
    print ('3. Check for harmonic constituents - required for harmonic & hydraulic stations')
    invalid_tide = tider.get_stations_with_invalid_harcons ()
    print_invalid_stations (invalid_tide, 'harmonic constituents', 'harmonic tide')
    print ('')
    invalid_current = currenter.get_stations_with_invalid_harcons ()
    print_invalid_stations (invalid_current, 'harmonic constituents', 'harmonic & hydraulic current')
    print ('')    

    ## Check datums for tides harmonic stations
    print ('4. Check for datums - required for harmonic tide stations')
    invalid_tide = tider.get_stations_with_invalid_datums ()
    print_invalid_stations (invalid_tide, 'datums', 'harmonic tide')
    print ('')

    ## Check major axis direction for current harmonic & hydraulic stations
    print ('5. Check for major axis direction - required for harmonic & hydraulic current stations')    
    invalid_current = currenter.get_stations_with_invalid_majordirs ()
    print_invalid_stations (invalid_current, 'major axis directions', 'harmonic & hydraulic current')
    print ('')

    ## Check reference stations for tides & current subordinate stations
    print ('6. Check for reference station - required for subordinate stations')
    invalid_tide = tider.get_stations_with_invalid_ref_info ()
    print_invalid_stations (invalid_tide, 'reference stations', 'subordinate tide')
    print ('')
    invalid_current = currenter.get_stations_with_invalid_ref_info ()
    print_invalid_stations (invalid_current, 'reference stations', 'subordinate current')
    print ('')

    ## Check prediction offsets for tides & current subordinate stations
    print ('6. Check for prediction offsets - required for subordinate stations')
    invalid_tide = tider.get_stations_with_invalid_offsets ()
    print_invalid_stations (invalid_tide, 'prediction offsets', 'subordinate tide')
    print ('')    
    invalid_current = currenter.get_stations_with_invalid_offsets ()
    print_invalid_stations (invalid_current, 'prediction offsets', 'subordinate current')
    print ('')

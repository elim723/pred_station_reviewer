#!C:/Users/elim.thompson/AppData/Local/Programs/Python/Python37/python

####
#### By Elim Thompson (01/28/2020)
####
#### The reviewer is a class to check API metadata for prediction stations. The child
#### classes - tide reviewer and current reviewer are inherited from the reviewer class.
####
#### For Harmonic stations ...
####  +-----------+---------------+--------------------+----------------------------------+
####  | Parameter | API           | Tide keys          | Current keys                     |
####  +-----------+---------------+--------------------+----------------------------------+
####  | time zone | stations.json | timezonecorr       | timezone_offset                  |
####  +-----------+---------------+--------------------+----------------------------------+
####  | harcon    | harcon.json   | number, amplitude  | binNbr, constNum, majorAmplitude |
####  |           |               | phase_local        | majorPhase, minorAmplitude       |
####  |           |               |                    | minorPhase, MajorMeanSpeed       |
####  |           |               |                    | MinorMeanSpeed                   |
####  +-----------+---------------+--------------------+----------------------------------+
####  | major dir |               | -                  | meanFloodDir                     |
####  +-----------+---------------+--------------------+----------------------------------+
####  | datums    | datums.json   | MSL                | n/a                              |
####  +-----------+---------------+--------------------+----------------------------------+
####
#### For Reference stations ...
####  +-----------+---------------+-----------------------+-------------------------------+
####  | Parameter | API           | Tide keys             | Current keys                  |
####  +-----------+---------------+-----------------------+-------------------------------+
####  | ref       | stations.json | refStationId          | refStationId, refStationBin   |
####  +-----------+---------------+-----------------------+-------------------------------+
####  | offsets   | stations.json | heightOffsetHighTide, | mfcTimeAdjMin, sbeTimeAdjMin  |
####  |           |               | heightOffsetLowTide,  | mecTimeAdjMin, sbfTimeAdjMin  |
####  |           |               | timeOffsetHighTide,   | mfcAmpAdj, mecAmpAdj          |
####  |           |               | timeOffsetLowTide,    |                               |
####  |           |               | heightAdjustedType    |                               |
####  +-----------+---------------+-----------------------+-------------------------------+
####
#### For Hydraulic stations (only current) ...
####  +-----------+---------------+--------------------+----------------------------------+
####  | Parameter | API           | Tide keys          | Current keys                     |
####  +-----------+---------------+--------------------+----------------------------------+
####  | time zone | stations.json | -                  | timezone_offset                  |
####  +-----------+---------------+--------------------+----------------------------------+
####  | harcon    | harcon.json   | -                  | binNbr, constNum, majorAmplitude |
####  |           |               |                    | majorPhase, minorAmplitude       |
####  |           |               |                    | minorPhase, MajorMeanSpeed       |
####  |           |               |                    | MinorMeanSpeed                   |
####  +-----------+---------------+--------------------+----------------------------------+
####  | major dir |               | -                  | meanFloodDir                     |
####  +-----------+---------------+--------------------+----------------------------------+
####  
##############################################################################################

import requests, numpy

###################################################
### Define Constants
###################################################
url = "https://tidesandcurrents.noaa.gov/mdapi/latest/webapi/stations"

nulls = ['', 'null', None]

tide_keys = {'timezone': 'timezonecorr',
             'harcon'  : ['number', 'amplitude', 'phase_local'],
             'datums'  : ['MSL'],
             'ref_info': 'refStationId',
             'offsets' : ['heightOffsetHighTide', 'heightOffsetLowTide', 'timeOffsetHighTide',
                          'timeOffsetLowTide'   , 'heightAdjustedType']}

current_keys = {'timezone': 'timezone_offset',
                'majordir': 'meanFloodDir',
                'harcon'  : ['binNbr', 'constNum', 'majorAmplitude', 'majorPhase', 'minorAmplitude',
                             'minorPhase', 'majorMeanSpeed', 'minorMeanSpeed'],
                'ref_info': ['refStationId', 'refStationBin'],
                'offsets' : ['mfcTimeAdjMin', 'sbeTimeAdjMin', 'mecTimeAdjMin',
                             'sbfTimeAdjMin', 'mfcAmpAdj'    , 'mecAmpAdj']}

###################################################
### Define Short Lambda Functions
###################################################
get_api = lambda url, sType, expand: "{0}.json?type={1}&expand={2}".format(url, sType, expand)
get_harcon = lambda url, station: "{0}/{1}/harcon.json".format(url, station)
get_datums = lambda url, station: "{0}/{1}/datums.json".format(url, station)

check_valid = lambda nulls, array: numpy.array ([v not in nulls for v in array])

###################################################
### Define Reviewer
###################################################
class reviewer (object):

    def __init__ (self, api):

        ''' Initialize a reviewer instance, and pull all data from the input API.

            input param
            -----------
            api (str): the API of the station meta data
        '''

        self.station_types = ['harmonic', 'subordinate', 'hydraulic']
        ## pull station metadata from API
        self._api   = api
        self._mdata = self._pull_metadata ()
        ## initialize private variables common to both tide and current prediction stations
        self.review_keys    = None # keys in metadata to be reviewed
        self._offset_apikey = None # tidepredictionoffsets / currentpredictionoffsets
        self._offsets       = None # offsets metadata for subordinate stations
        self._ref_info      = None # reference stations for subordinate stations

    @property
    def harmonic_stations (self):
        ''' Getter for sorted harmonic stations '''
        return numpy.array (sorted (self._mdata['harmonic']))

    @property
    def subordinate_stations (self):
        ''' Getter for sorted subordinate stations '''
        return numpy.array (sorted (self._mdata['subordinate']))

    @property
    def hydraulic_stations (self):
        ''' Getter for sorted hydraulic stations '''
        return numpy.array (sorted (self._mdata['hydraulic']))

    @property
    def offsets (self):
        ''' Getter for offset data for subordinate stations '''
        return self._offsets

    @property
    def ref_info (self):
        ''' Getter for reference stations for subordinate stations '''
        return self._ref_info

    def _pull_metadata (self):

        ''' Private function to pull station metadata and organize data in types.
            Only 3 types are included: harmonic, subordinate, and hydraulic. Tide
            prediction stations are either harmonic or subordinate. For current,
            it can be harmonic, subordinate, hydraulic, or weak. No predictions
            are generated for weak current stations, so they are excluded from
            metadata reviewal. A master dictionary is returned with all meta-
            data with key of station and value of the corresponding metadata.

            return param
            ------------
            mdict (dict): master dictionary containing all station metadata
        '''

        ## pull the 'stations' metadata from API
        response = requests.get (self._api)
        mdata = response.json ()
        mdata = mdata['stations']
        ## group stations by type: harmonic, subodinate, and hydraulic.
        ## No hydraulic tide prediction stations is expected.
        mdict = {stype:{} for stype in self.station_types}
        for jdata in mdata:
            ## tide station is 9 digit e.g. 9450623
            ## current station has both ID and Bin e.g. ACT1616_1
            station = jdata["id"]
            if "currbin" in jdata.keys():
                station += '_' + str (jdata["currbin"]) 
            ## only study type R / H (harmonic), S (subordinate), and Y (hydraulic)
            stype = 'harmonic'    if jdata["type"] in ['R', 'H'] else \
                    'subordinate' if jdata["type"] == 'S' else \
                    'hydraulic'   if jdata["type"] == 'Y' else None
            ## put this station metadata to master dictionary
            if stype is not None: mdict[stype][station] = jdata

        return mdict

    def _pull_harcon (self, station):
        
        ''' Private function to pull harmonic constituent data from station API. This
            function access the station/harcon.json API for the input station. The data
            is then re-organized into a dictionary. Its key are the same as the ones in
            self.review_keys['harcon']. Each key is an array with the same size as # of
            constituents for that station. For current, it is the same except that only
            the constituents for the input station bin are included.

            input param
            -----------
            station (string): station to pull harcon. For tide, it is 9 digit string e.g.
                              '8517756'. Current ones must have bin info e.g. 'ACT1616_1'.
            
            return param
            ------------
            harcon_dict (dict): Constituent values for the input station. Keys are
                                the same as self.review_keys['harcon'], values must
                                have length the same as # constituents available.
        '''

        ## Keys where data is extracted and reviewed
        harcon_keys = self.review_keys['harcon']

        ## pull harcon data by looking into http://URL/{station}/harcon.json.
        #  for current, {station} on API does not include bin info.
        station_id = station
        if '_' in station: station_id, station_bin = station.split ('_')
        #  extract all constituent data
        response = requests.get (get_harcon (url, station_id))
        mdata = response.json ()
        harcon = mdata['HarmonicConstituents']
        #  if none is available, return an empty dictionary
        if len (harcon) == 0: return {key:[] for key in harcon_keys}

        ## organize harcon data based on harcon_keys
        harcon_data = [[consti[key] for key in harcon_keys] for consti in harcon]
        harcon_data = numpy.array (harcon_data).T
        harcon_dict = {key:array for key, array in zip (harcon_keys, harcon_data)}
        #  for current, select the constituent values for this specific bin.
        if 'binNbr' in harcon_dict: 
            is_this_bin = harcon_dict['binNbr'] == int (station_bin)
            harcon_dict = {key:array[is_this_bin] for key, array in harcon_dict.items()}
        return harcon_dict

    def _print_info (self, n_invalid, info_name, stype, do_print=False):

        ''' Private function to show # invalid stations on console. Print when asked.

            input param
            -----------
            n_invalid (int)   : number of stations with invalid information
            info_name (string): name of the invalid information
            stype     (string): station type e.g. tide harmonic / current hydraulic / etc
            do_print  (bool)  : If True, print the message.
        '''
        ## Print only when asked
        if not do_print: return
        ## Print the message.
        message = "{0} {1} stations have invalid {2}.".format(n_invalid, stype, info_name)
        print (message)

    def _collect_timezones (self, stype):
        ''' Private function to extract time zones. Only applied if input station type is either
            harmonic or hydraulic. Time zone value is pulled directly from the main station.json
            API. No need to access the individual station API http://URL/{station}.json.

            input param
            -----------
            stype (str): station type - either 'harmonic' or 'hydraulic'
            
            return param
            ------------
            timezones (array): time zones of the input station type.
        '''
        ## Only allow harmonic / hydraulic station type.
        if not stype in ['harmonic', 'hydraulic']:
            raise IOError ('Only harmonic / hydraulic are allowed for time zones. No {0}.'.format(stype))

        ## Loop through the main metadata and extract time zones for the specific station type.
        stations = self.harmonic_stations if stype == 'harmonic' else self.hydraulic_stations
        timezones = [self._mdata[stype][station][self.review_keys['timezone']] for station in stations]
        return numpy.array (timezones)

    def _collect_harcon (self, stype):
        ''' Private function to extract harmonic constituent values. Only applied if input station
            type is either harmonic or hydraulic. Constituent values are pulled from the individual
            station API http://URL/{station}.json. This function loops through each station of the
            requested type to collect all harmonic constituents for all stations.

            The output dictionary has keys equal to the station ID (9-digit string for tide and with
            bin info for current). The values are sub-dictionary from self._pull_harcon function.

            input param
            -----------
            stype (str): station type - either 'harmonic' or 'hydraulic'
            
            return param
            ------------
            harcon_dict (dict): harmonic constituent values of the input station type.
        '''
        ## Only allow harmonic / hydraulic station type.
        if not stype in ['harmonic', 'hydraulic']:
            raise IOError ('Only harmonic / hydraulic are allowed for harmonic constituents. No {0}.'.format(stype))

        ## Loop through the main metadata and extract time zones for the specific station type.
        stations = self.harmonic_stations if stype == 'harmonic' else self.hydraulic_stations
        return {station:self._pull_harcon (station) for station in stations}

    def _collect_ref_info (self, is_tide=False):
        ''' Private function to map subordinate stations with their reference stations. Only applied
            to subordinate stations. Their reference station are pulled from the main station.json
            API, specifically the self._offset_apikey and self.review_keys['ref_info']. No need to
            access the individual station API http://URL/{station}.json.
            
            This function loops through each subordinate station to collect all their reference
            station. For tide, stations are named as a 9-digit string, whereas current stations have
            bin information in their station ID. These sub-ref pairs are re-organized as a dictionary
            where the keys are reference stations, and each value is an array of all subordinate
            stations that reference to the key reference station.

            input param
            -----------
            is_tide (bool): If True, station ID are expected to be 9-digit strings. If False, station
                            ID has two pieces - ID + Bin number.
            
            return param
            ------------
            ref_dict (dict): Keys are the harmonic stations that are referenced by at least one
                             subordinate stations. Each key has an array value listing all
                             subordinate station that reference to the key reference station.
        '''

        ## Keys where data is extracted and reviewed
        ref_keys = self.review_keys['ref_info']

        ## Loop through each subordinate stations. Put all reference stations into a dictionary.
        ## Dictionary is organized by the reference station each of which has a list of its
        ## subordinate stations.
        ref_dict = {}
        for sub in self.subordinate_stations:
            # extract the prediction offsets metadata for this subordinate station
            ref_mdata = self._mdata['subordinate'][sub][self._offset_apikey]
            # handle the difference between tide station ID and current ID
            ref = ref_mdata[self.review_keys['ref_info']] if is_tide else \
                  ref_mdata[self.review_keys['ref_info'][0]] + '_' + str (ref_mdata[self.review_keys['ref_info'][1]])
            # if this ref station already exists in ref_dict, append the sub station to its list.
            # otherwise, create a new key for this ref station.
            if ref in ref_dict.keys():
                ref_dict[ref] += [sub]
            else:
                ref_dict[ref] = [sub]
        ## Return the same dictionary but with all values numpified.
        return {ref:numpy.array (sub) for ref, sub in ref_dict.items()}

    def _collect_offsets (self):
        ''' Private function to extract prediction offset values. Only applied to subordinate stations.
            Their offset values are from the main station.json API, specifically the self._offset_apikey
            and self.review_keys['offsets']. No need to access the individual station API /{station}.json.
            
            This function loops through each subordinate station to collect all their offset values.
            For each key in self.review_keys['offsets'], we expect 1 value from each subordinate stations.
            Therefore, these offsets data are re-organized as a dictionary where the keys are the same as
            self.review_keys['offsets'], and each value is an array (size = number of subordinate stations)
            of the corresponding offset types.
            
            return param
            ------------
            offset_dict (dict): Keys are the same as self.review_keys['offsets']. Each value is an array
                                (size = # subordinate stations) of the corresponding offset type.
        '''

        ## Keys where data is extracted and reviewed
        offset_keys = self.review_keys['offsets']
        ## Pull all subordinate metadata
        sdata = self._mdata['subordinate']
        ## Loop through each subordinate station and pull all the offset keys for reviewal
        offset_data = [[sdata[sub][self._offset_apikey][key] for key in offset_keys]
                       for sub in self.subordinate_stations]
        offset_data = numpy.array (offset_data).T
        ## Re-format offset data to a dictionary
        return {key:array for key, array in zip (offset_keys, offset_data)}

    def _check_data_validity (self, array, info_name='', stype='', do_print=False):
    
        ''' Private function to check for any invalid data in the input array. An element in the
            array is invalid if it is equal to one of the null values i.e. '', 'null', or None. 

            input param
            -----------
            array   (array): numpy array to be checked
            info_name (str): name of the information being checked
            stype     (str): station type e.g. harmonic tidal / subordinate current / etc.
            do_print (bool): If True, print out how many stations are invalid.

            return param
            ------------
            valid (array): a boolean array indicating which element is valid in the array.
        '''
        ## check which elements are valid and which are not
        valid = check_valid (nulls, array)
        ## print out how man invalid data is found from this array
        self._print_info (len (valid[~valid]), info_name, stype, do_print=do_print)
        return valid

    def _get_stations_with_invalid_offsets (self, stype, do_print=False):
        ''' Private function to return an array of subordinate stations with invalid offset
            values. For each subordinate station, we make sure all the offsets in
            self.review_keys['offsets'] are available. For 'heightAdjustedType', its value
            can either be 'R' or 'F' - all other values are considered as invalid. These
            offsets values were already pulled from the main API - no need to access
            individual station APIs. If do_print is True, this function will let you know
            how many subordinate stations with invalid offsets data are found. This function
            is only applicable to subordinate stations. This function is shared between the
            tide & current stations.

            input param
            -----------
            stype     (str): station type for printing messages
            do_print (bool): If True, print # subordinate stations with invalid offsets data.

            return param
            ------------
            station_array (array): an array of subordinate stations with invalid offsets data.
        '''
        ## Keys where data is extracted and reviewed
        offsets_keys = self.review_keys['offsets']

        ## If self._offsets is not set, collect the offsets data now.
        if self.offsets is None: self._offsets = self._collect_offsets ()
        ## Check each of the sum key. 'HeightAdjustedType' can either be 'R' or 'F', while the
        ## rest are just numbers. A station has invalid offsets if any one of those keys have
        ## invalid value. The valid variable is an array of size the same as # subordinate
        ## stations and keeps track of the number of valid offsets at each station. Any number
        ## less than len (self.review_keys['offsets']) means that at least one of the offset
        ## data for this subordinate station is invalid.
        valid = numpy.zeros (len (self.subordinate_stations))
        for key in offsets_keys:
            are_valid = [v in ['R', 'F'] for v in self.offsets[key]] if key == 'heightAdjustedType' else \
                        self._check_data_validity (self.offsets[key], do_print=False)       
            valid += numpy.array (are_valid).astype (int)
        invalid_offsets = valid < len (offsets_keys)
        ## Let user know of any subordinate stations with invalid offsets
        self._print_info (len (invalid_offsets[invalid_offsets]), 'offsets', stype=stype, do_print=do_print)
        return self.subordinate_stations[invalid_offsets]

###################################################
### Define Tide Reviewer
###################################################
class tide_reviewer (reviewer):
    
    def __init__ (self):
        ''' Initialize a tide reviewer inherited from reviewer class. Three pieces of data
            are unique to tide prediction stations:
                * (harmonic stations) datums
                * (harmonic stations) time zones
                * (harmonic stations) harmonic constituents
            No hydraulic tide stations are reviewed.
        '''
        ## Inherit from reviewer class.
        super().__init__ (get_api (url, "tidepredictions", "tidepredoffsets"))
        ## initialize private variables specific to tide stations
        self.review_keys    = tide_keys         # metadata to be reviewed
        self._offset_apikey = 'tidepredoffsets' # prediction offsets keys on API
        self._datums        = None              # unique metadata for tide
        self._timezones     = None              # unique metadata for tide
        self._harcons       = None              # unique metadata for tide

    @property
    def datums (self):
        ''' Getter for datums data for harmonic tide stations '''
        return self._datums

    @property
    def timezones (self):
        ''' Getter for time zone data for harmonic tide stations '''
        return self._timezones

    @property
    def harcons (self):
        ''' Getter for harmonic constituent data for harmonic tide stations '''
        return self._harcons

    @ property
    def bad_ref_stations (self):
        ''' Getter for all bad reference stations '''
        return numpy.array (list (self.get_stations_with_invalid_timezones (do_print=False)) + \
                            list (self.get_stations_with_invalid_harcons (do_print=False)) + \
                            list (self.get_stations_with_invalid_datums (do_print=False)))

    def _pull_datums (self, station):
        ''' Private function to pull datums data from station API. This function access the
            station/datums.json API for the input station. All datum values are first pulled.
            These data is then re-organized into an array in the order of datum names in
            self.review_keys['datums']. This function is only applicable to harmonic stations
            since no datums is required for subordinate tide stations.

            input param
            -----------
            station (int): harmonic tide station to pull harcon. 
            
            return param
            ------------
            datums_array (array): Datum values in the order of datum names defined in
                                  self.review_keys['datums']
        '''
        ## Pull datums data from this station API 
        response = requests.get (get_datums (url, station))
        mdata = response.json ()
        datums = mdata['datums']
        ## Re-format the datums data
        datums_dict = {datum['name']:datum['value'] for datum in datums}
        ## Return those to be reviewed as an array
        return [datums_dict[key] if key in datums_dict.keys() else None
                for key in self.review_keys['datums']]

    def _collect_datums (self):
        ''' Private function to extract datum values. Only applied if input station type is harmonic.
            Datum values are pulled from the individual station API http://URL/{station}.json. This
            function loops through each harmonic station to collect the datum values to be reviewed.

            The output dictionary has keys equal to the station ID (9-digit string for tide). Each key
            has a value of type array storing the datum values for that station ID in the ordering in
            self.review_keys['datums'].
            
            return param
            ------------
            datum_dict (dict): Each key is a station ID. Each value is an array of datums for that
                               station to be reviewed.
        '''
        ## Collect datums for reviewal from all harmonic stations
        return {station:self._pull_datums (station) for station in self.harmonic_stations}

    def get_stations_with_invalid_timezones (self, do_print=False):
        ''' Public function to return an array of harmonic stations with invalid time zones.
            A time zone is considered as invalid if its value is either '', 'null', or None.
            If do_print is True, this function will let you know how many harmonic stations
            with invalid time zones are found. This function is only applicable to harmonic
            stations - subordinate stations do not require time zone to be valid.

            input param
            -----------
            do_print (bool): If True, print # harmonic stations with invalid time zone.

            return param
            ------------
            station_array (array): an array of harmonic stations with invalid time zones.
        '''

        ## If self._timezones is not set, collect all time zone data now.
        if self.timezones is None: self._timezones = self._collect_timezones ('harmonic')
        ## Check validity of time zones for all harmonic stations
        valid = self._check_data_validity (self.timezones, info_name='time zones',
                                           stype='harmonic tidal', do_print=do_print)
        ## return harmonic stations with invalid time zones
        return self.harmonic_stations[~valid]
    
    def get_stations_with_invalid_harcons (self, do_print=False):
        ''' Public function to return an array of harmonic stations with invalid constituent
            data. Constituent data comes as sets. If invalid, the set has a length of 0.0 i.e.
            no constituents are available on API. If do_print is True, this function will let
            you know how many harmonic stations with invalid constituents are found. This
            function is only applicable to harmonic stations - subordinate stations do not
            require constituent data to be valid.

            input param
            -----------
            do_print (bool): If True, print # harmonic stations with invalid constituent data.

            return param
            ------------
            station_array (array): an array of harmonic stations with invalid constituent data.
        '''
        ## If self._harcon is not set, collect all constituent data now.
        if self.harcons is None: self._harcons = self._collect_harcon ('harmonic')
        ## Check validity based on the number of constituents found for each station.
        valid = numpy.array ([len (harcon[self.review_keys['harcon'][0]]) > 0
                              for station, harcon in self.harcons.items()])
        ## Print out a message if asked
        self._print_info (len (valid[~valid]), 'harmonic constituents',
                          'harmonic tidal', do_print=do_print)
        ## return harmonic stations with no constituents on API
        return self.harmonic_stations[~valid]

    def get_stations_with_invalid_datums (self, do_print=False):
        ''' Public function to return an array of harmonic stations with invalid datum values.
            A datum value is considered as invalid if it is either '', 'null', or None. If
            do_print is True, this function will let you know how many harmonic stations with
            invalid datums are found. This function is only applicable to harmonic stations -
            subordinate stations do not require constituent data to be valid. A minimal datum
            required for prediction calculation is MSL.

            input param
            -----------
            do_print (bool): If True, print # harmonic stations with invalid datums data.

            return param
            ------------
            station_array (array): an array of harmonic stations with invalid datums data.
        '''        
        ## If self._datums is not set, collect all datums data now.
        if self.datums is None: self._datums = self._collect_datums ()
        ## Check validity of datums for all harmonic stations
        valid = self._check_data_validity (self.datums, info_name='datums',
                                           stype='harmonic tidal', do_print=do_print)        
        ## return harmonic stations with invalid datums
        return self.harmonic_stations[~valid]

    def get_stations_with_invalid_ref_info (self, do_print=False):
        ''' Public function to return a dictionary of invalid reference stations and their
            subordinate stations. For each subordinate station, we check for the validity
            of its reference station:
                1. Its ref station must be a harmonic station
                2. This harmonic station has valid time zone, harcon, and datums.
            If do_print is True, this function will let you know how many subordinate stations
            with invalid ref info are found. This function is only applicable to subordinate
            stations.

            input param
            -----------
            do_print (bool): If True, print # subordinate stations with invalid ref data.

            return param
            ------------
            ref-sub (dict): a dictionary where keys are bad reference station, and their values
                            are ararys of the corresponding subordinate stations            
        '''
        ## If self._ref_info is not set, build a ref-sub map now.
        if self.ref_info is None: self._ref_info = self._collect_ref_info (is_tide=True)
        ## For each ref station, check if ..
        ##   1. it is a harmonic station
        ##   2. it has valid time zone, harcon, and datums.
        valid_ref_station = numpy.array ([ref in self.harmonic_stations and ref not in self.bad_ref_stations
                                          for ref in self.ref_info.keys()])
        ## Collect subordinate stations that reference to bad ref stations
        invalid_ref_stations = numpy.array (list (self.ref_info.keys()))[~valid_ref_station]
        n_subordinate_with_invalid_ref = len ([sub for ref in invalid_ref_stations for sub in self.ref_info[ref]])
        ## Print message if asked                                      
        self._print_info (n_subordinate_with_invalid_ref, 'ref stations', 'subordinate tidal', do_print=do_print)
        ## Return a map of bad ref (key) - sub array (value)
        return {ref:self.ref_info[ref] for ref in invalid_ref_stations}

    def get_stations_with_invalid_offsets (self, do_print=False):
        ''' Public function to return an array of subordinate stations with invalid offset
            values. For each subordinate station, we make sure all the offsets in
            self.review_keys['offsets'] are available. For 'heightAdjustedType', its value
            can either be 'R' or 'F' - all other values are considered as invalid. These
            offsets values were already pulled from the main API - no need to access
            individual station APIs. If do_print is True, this function will let you know
            how many subordinate stations with invalid offsets data are found. This function
            is only applicable to subordinate stations.

            input param
            -----------
            do_print (bool): If True, print # subordinate stations with invalid offsets data.

            return param
            ------------
            station_array (array): an array of subordinate stations with invalid offsets data.
        '''
        return self._get_stations_with_invalid_offsets ('subordinate tidal', do_print=do_print)

###################################################
### Define Current Reviewer
###################################################
class current_reviewer (reviewer):
    
    def __init__ (self):
        ''' Initialize a current reviewer inherited from reviewer class. Three pieces of data
            are unique to current prediction stations:
                 * (harmonic stations) time zones
                 * (harmonic stations) harmonic constituents
                 * (harmonic stations) major axis direction
                 * (hydraulic stations) time zones
                 * (hydraulic stations) harmonic constituents
                 * (hydraulic stations) major axis direction
        '''
        ## Inherit from reviewer class.
        super().__init__ (get_api (url, "currentpredictions", "currentpredictionoffsets"))
        ## initialize private variables specific to tide stations
        self.review_keys    = current_keys               # metadata to be reviewed
        self._offset_apikey = 'currentpredictionoffsets' # prediction offsets keys on API
        self._harmonic_timezones = None
        self._harmonic_harcons   = None
        self._harmonic_majordirs = None
        self._hydraulic_timezones = None
        self._hydraulic_harcons   = None
        self._hydraulic_majordirs = None

    @property
    def harmonic_timezones (self):
        ''' Getter for time zone data for harmonic current stations '''
        return self._harmonic_timezones

    @property
    def hydraulic_timezones (self):
        ''' Getter for time zone data for hydraulic current stations '''
        return self._hydraulic_timezones

    @property
    def harmonic_harcons (self):
        ''' Getter for harmonic constituent data for harmonic current stations '''
        return self._harmonic_harcons

    @property
    def hydraulic_harcons (self):
        ''' Getter for harmonic constituent data for hydraulic current stations '''
        return self._hydraulic_harcons

    @property
    def harmonic_majordirs (self):
        ''' Getter for major axis direction data for harmonic current stations '''
        return self._harmonic_majordirs

    @property
    def hydraulic_majordirs (self):
        ''' Getter for major axis direction data for hydraulic current stations '''
        return self._hydraulic_majordirs

    @ property
    def bad_ref_stations (self):
        ''' Getter for all bad reference stations '''
        return numpy.array (list (self.get_stations_with_invalid_timezones (do_print=False)) + \
                            list (self.get_stations_with_invalid_harcons (do_print=False)) + \
                            list (self.get_stations_with_invalid_majordirs (do_print=False)))

    def _collect_majordirs (self, stype):
        ''' Private function to gather major axis directions. Only applied if input station type is either
            harmonic or hydraulic. Major axis directions are pulled directly from the main station API
            station.json. This function loops through each harmonic / hydraulic station to collect the
            major axis directions.

            Because one harmonic / hydraulic station has one direction value, the output is an array
            with the size as the # stations of requested type.
            
            input param
            -----------
            stype (str): station type - either 'harmonic' or 'hydraulic'

            return param
            ------------
            majordirs (array): Major axis directions for stations of requested type.
        '''
        
        ## Make sure input stype is valid
        if not stype in ['harmonic', 'hydraulic']:
            raise IOError ('{0} does not require valid time zone; Only harmonic / hydraulic.'.format(stype))

        ## Collect and return all major axis direction for the specific station type
        stations = self.harmonic_stations if stype == 'harmonic' else self.hydraulic_stations
        majordirs = [self._mdata[stype][station][self._offset_apikey][self.review_keys['majordir']]
                     for station in stations]
        return numpy.array (majordirs)

    def get_stations_with_invalid_timezones (self, do_print=False):
        ''' Public function to return an array of harmonic & hydraulic stations with invalid
            time zones. A time zone is considered as invalid if its value is either '', 'null',
            or None. If do_print is True, this function will let you know how many harmonic /
            hydraulic stations with invalid time zones are found. This function is only
            applicable to harmonic & hydraulic stations - subordinate stations do not require
            time zone to be valid.

            input param
            -----------
            do_print (bool): If True, print the numbers of harmonic & hydraulic stations with
                             invalid time zone. One message per station type.

            return param
            ------------
            station_array (array): an array of harmonic + hydraulic stations with invalid time
                                   zones. One list with both station types is returned.
        '''

        ## initialize a holder to store harmonic & hydraulic stations with invalid time zones.
        invalid_stations = []
        ## loop through each station type
        for stype in ['harmonic', 'hydraulic']:
            # get the time zone values for this station type. If it is not yet exist, collect
            # the time zone data now.
            this_attr = '_' + stype + '_timezones'
            if getattr (self, this_attr, None) is None:
                setattr (self, this_attr, self._collect_timezones (stype))   
            timezones = getattr (self, this_attr)
            # check if time zone array is valid
            valid = self._check_data_validity (timezones, info_name='time zones',
                                               stype=stype + ' current', do_print=do_print)
            # append those with invalid time zones to invalid_stations holder
            stations = self.harmonic_stations if stype == 'harmonic' else self.hydraulic_stations
            invalid_stations += list (stations[~valid])
        return numpy.array (invalid_stations)

    def get_stations_with_invalid_majordirs (self, do_print=False):
        ''' Public function to return an array of harmonic & hydraulic stations with invalid
            major axis direction. A direction is considered as invalid if its value is either
            '', 'null', or None. If do_print is True, this function will let you know how many
            harmonic / hydraulic stations with invalid direction are found. This function is
            only applicable to harmonic & hydraulic stations - subordinate stations do not
            require major axis direction to be valid.

            input param
            -----------
            do_print (bool): If True, print the numbers of harmonic & hydraulic stations with
                             invalid major axis direction. One message per station type.

            return param
            ------------
            station_array (array): an array of harmonic + hydraulic stations with invalid major
                                   axis direction. One list with both station types is returned.
        '''
        ## initialize a holder to store harmonic & hydraulic stations
        ## with invalid major axis directions.
        invalid_stations = []
        ## loop through each station type
        for stype in ['harmonic', 'hydraulic']:
            # get the major axis direction for this station type. If it is not yet exist, collect
            # the major axis direction now.
            this_attr = '_' + stype + '_majordirs'
            if getattr (self, this_attr, None) is None:
                setattr (self, this_attr, self._collect_majordirs (stype))
            majordirs = getattr (self, this_attr)
            # check if major axis direction array is valid
            valid = self._check_data_validity (majordirs, info_name='major axis direction',
                                               stype=stype + ' current', do_print=do_print)
            # append those with invalid major axis direction to invalid_stations holder
            stations = self.harmonic_stations if stype == 'harmonic' else self.hydraulic_stations
            invalid_stations += list (stations[~valid])
        return numpy.array (invalid_stations)
    
    def get_stations_with_invalid_harcons (self, do_print=False):
        ''' Public function to return an array of harmonic & hydraulic stations with invalid
            constituent data. Constituent data comes as sets. If invalid, the set has a length
            of 0.0 i.e. no constituents are available on API. If do_print is True, this function
            will let you know how many harmonic / hydraulic stations with invalid constituents
            are found. This function is only applicable to harmonic & hydraulic stations.
            Subordinate stations do not require constituent data to be valid.

            input param
            -----------
            do_print (bool): If True, print the number of harmonic / hydraulic stations with
                             invalid constituent data. One message per station type.

            return param
            ------------
            station_array (array): an array of harmonic & hydraulic stations with invalid
                                   constituent data. One array for both station types.
        '''        
        ## initialize a holder to store harmonic & hydraulic stations with no constituents.
        invalid_stations = []
        ## loop through each station type
        for stype in ['harmonic', 'hydraulic']:
            # get the constituent data for this station type. If it is not yet exist, collect
            # the constituent data now.
            this_attr = '_' + stype + '_harcons'
            if getattr (self, this_attr, None) is None:
                setattr (self, this_attr, self._collect_harcon (stype))
            harcons = getattr (self, this_attr)
            ## Check validity based on the number of constituents found for each station.
            valid = numpy.array ([len (harcon[self.review_keys['harcon'][0]]) > 0
                                  for station, harcon in harcons.items()])            
            ## Print out a message if asked
            self._print_info (len (valid[~valid]), 'harmonic constituents',
                              stype + ' current', do_print=do_print)
            # append those with no constituents to invalid_stations holder
            stations = self.harmonic_stations if stype == 'harmonic' else self.hydraulic_stations
            invalid_stations += list (stations[~valid])
        return numpy.array (invalid_stations)

    def get_stations_with_invalid_ref_info (self, do_print=False):
        ''' Public function to return a dictionary of invalid reference stations and their
            subordinate stations. For each subordinate station, we check for the validity
            of its reference station:
                1. Its ref station must be a harmonic / hydraulic station
                2. This harmonic / hydraulic station has valid time zone, harcon, and datums.
            If do_print is True, this function will let you know how many subordinate stations
            with invalid ref info are found. This function is only applicable to subordinate
            stations.

            input param
            -----------
            do_print (bool): If True, print # subordinate stations with invalid ref data.

            return param
            ------------
            ref-sub (dict): a dictionary where keys are bad reference station, and their values
                            are ararys of the corresponding subordinate stations
        '''
        ## If self._ref_info is not set, build a ref-sub map now.
        if self.ref_info is None: self._ref_info = self._collect_ref_info (is_tide=False)
        ## For each ref station, check if ..
        ##   1. it is a harmonic / hydraulic station
        ##   2. it has valid time zone, harcon, and datums.
        valid_ref_station = numpy.array ([(ref in self.harmonic_stations or ref in self.hydraulic_stations)
                                          and ref not in self.bad_ref_stations for ref in self.ref_info.keys()])
        ## Count how many subordinate stations with bad ref stations
        invalid_ref_stations = numpy.array (list (self.ref_info.keys()))[~valid_ref_station]
        n_subordinate_with_invalid_ref = len ([sub for ref in invalid_ref_stations for sub in self.ref_info[ref]])
        ## Print message if asked 
        self._print_info (n_subordinate_with_invalid_ref, 'ref stations',
                          'subordinate current', do_print=do_print)
        ## Return a map of bad ref (key) - sub array (value)
        return {ref:self.ref_info[ref] for ref in invalid_ref_stations}

    def get_stations_with_invalid_offsets (self, do_print=False):
        ''' Public function to return an array of subordinate stations with invalid offset
            values. For each subordinate station, we make sure all the offsets in
            self.review_keys['offsets'] are available. These offsets values were already
            pulled from the main API - no need to access individual station APIs. If do_print
            is True, this function will let you know how many subordinate stations with
            invalid offsets data are found. This function is only applicable to subordinate
            stations.

            input param
            -----------
            do_print (bool): If True, print # subordinate stations with invalid offsets data.

            return param
            ------------
            station_array (array): an array of subordinate stations with invalid offsets data.
        '''
        return self._get_stations_with_invalid_offsets ('subordinate current', do_print=do_print)
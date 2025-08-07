import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

class NS_Solver:
    """
    A class to solve the quadratic tilt and deflection of a bridge deck using ascending and descending displacement data.
    This class initializes with data containing ascending and descending displacement information,
    sets up polynomial functions for displacement, and provides methods to calculate tilt and deflection.   
    
    Attributes
    ----------
    data : dict
        A dictionary containing ascending and descending displacement data.
    polyfunction : dict
        A dictionary containing polynomial functions for ascending and descending displacement.

    Methods
    -------
    setup() -> None:
        Sets up the polynomial functions for displacement.
    quadratic_tilt(keyword: str) -> float:
        Calculates the quadratic tilt of the bridge deck.
    quadratic_deflection(keyword: str) -> float:
        Calculates the quadratic deflection of the bridge deck.
    _quadratic_x(orbit: str) -> np.ndarray:
        Generates a linear space of x values for the specified orbit.
    _quadratic_y(orbit: str) -> np.ndarray:
        Evaluates the polynomial function for the specified orbit.
    analytical_curve(orbit: str) -> np.ndarray:
        Computes the analytical curve for the specified orbit based on the polynomial fit.
    

    """
    def __init__(self, data):
        """
        Initialize the NS_Solver with the provided data.
        
        Args:
            data (dict): Data required for the solver.
        """
        self.data = data
        self.setup()
    
    def setup(self):
        """
        Set up the solver with the necessary configurations.
        This method can be overridden by subclasses to provide specific setup logic.
        """
        self.polyfunction = dict(
            ascending = np.poly1d(
                np.polyfit(
                    x = self.data['ascending']['ndist'],
                    y = self.data['ascending']['disp'],
                    deg = 2
                )
            ),
            descending = np.poly1d(
                np.polyfit(
                    x = self.data['descending']['ndist'],
                    y = self.data['descending']['disp'],
                    deg = 2
                )
            )
        )

    def quadratic_tilt(self , keyword:str) -> float:
        """ Calculate the quadratic tilt of the bridge deck.

        This method computes the tilt based on the polynomial fit of the displacement data.
        
        Args:
            keyword (str): The keyword to identify the dataset ('ascending' or 'descending').
        Returns:
            float: The maximum tilt ratio of the bridge deck.
        """
        xleft = self.data[keyword]['ndist'][0]
        xright = self.data[keyword]['ndist'][-1]
        
        return float(np.abs(self.polyfunction[keyword](xright) - self.polyfunction[keyword](xleft)) / self.data['deck']['deck_length'][0])

    def quadratic_deflection(self, keyword:str) -> float:
        """ Calculate the quadratic deflection of the bridge deck.
        
        This method computes the deflection based on the polynomial fit of the displacement data.

        Args:
            keyword (str): The keyword to identify the dataset ('ascending' or 'descending').
        Returns:
            float: The maximum deflection ratio of the bridge deck.
        """
        xleft, xright = self.data[keyword]['ndist'][0], self.data[keyword]['ndist'][-1]
        yleft, yright = self.polyfunction[keyword](xright), self.polyfunction[keyword](xleft)

        slope = (yright - yleft) / (xright - xleft)
        intercept = yleft - slope * xleft
        
        # xrange = np.linspace(xleft, xright, self.data[keyword]['ndist'].shape[0])
        xrange = self._quadratic_x(keyword)
        yinterp = self.polyfunction[keyword](xrange)
        
        deflection = yinterp - (slope * xrange + intercept) 
        return float(np.abs(deflection).max() / self.data['deck']['deck_length'][0])

    def _quadratic_x(self, orbit:str) -> np.ndarray:
        """ Generate a linear space of x values for the specified orbit.

        This method creates a linear space of x values based on the normalized distance of the specified orbit.
        
        Arguments
        ---------
        orbit : str
            The orbit type ('ascending' or 'descending').
        
        Returns
        -------
        np.ndarray: A linear space of x values for the specified orbit.
        """
        ndist =  self.data[orbit]['ndist']
        return np.linspace(ndist.min(), ndist.max(), 50)
    
    def _quadratic_y(self, orbit:str) -> np.ndarray:
        """ Evaluate the polynomial function for the specified orbit.
        
        This method evaluates the polynomial function for the specified orbit using the normalized distance.
        
        Arguments
        ---------
        orbit : str
            The orbit type ('ascending' or 'descending').
        
        Returns
        -------
        np.ndarray: The evaluated polynomial function values for the specified orbit.
        """
        ndist =  self._quadratic_x(orbit)
        return self.polyfunction[orbit](ndist)
    
    def analytical_curve(self, orbit:str) -> np.ndarray:
        """ Compute the analytical curve for the specified orbit based on the polynomial fit.

        This method generates a linear space of x values and evaluates the polynomial function for the specified orbit.

        Arguments
        ---------
            orbit (str): The orbit type ('ascending' or 'descending').
        Returns
        -------
            np.ndarray: The evaluated polynomial function values for the specified orbit.
        
        Methods
        -------
        fit_beam_displacement(x: np.ndarray, y: np.ndarray, L0: float) -> np.ndarray:
            Fits the beam displacement for a single span.
        
        one_span_beam_displacement(x: np.ndarray, L: float, C0: float, A: float, B: float) -> np.ndarray:
            Calculates the single span beam displacement.
        fit_beam_displacement_two_spans(x: np.ndarray, y: np.ndarray, L0: float, p0: np.ndarray = None, maxfev: int = 500) -> np.ndarray:
            Fits the beam displacement for two spans.
        two_span_beam_displacement(x: np.ndarray, L: float, C0: float, A: float, B: float, C: float) -> np.ndarray:
            Calculates the two span beam displacement.
        """
        
        deck_length = self.data['deck']['deck_length'][0]
        span_count = self.data['deck']['span_count'][0]
        xrange = self._quadratic_x(orbit)
        
        def fit_beam_displacement(x:np.ndarray, y:np.array, L0:float) -> np.ndarray:
            """ Method fits the single span beam displacement
            
            Arguments
            -----------
            x : np.array
                distance from the selected point along the bridge
            y : np.array
                displacement values
            L0 : float
                length of the deck
            
            Returns
            --------
            opt_params: optimized parameters for the single span beam displacement
            """
            opt_params, _ = curve_fit(lambda x, C0, A, B: one_span_beam_displacement(x, L0, C0, A, B), x, y)
            return opt_params
        
        def one_span_beam_displacement(x:np.ndarray, L:float, C0:float, A:float, B:float) -> np.ndarray:
            """method calculates and returns the single span beam displacement
        
            Arguments
            -----------
            x : np.array
                distance from the selected point along the bridge
            L : float
                length of the deck
            C0, A, B : float
                interpolation variables
            
            Returns
            --------
            np.ndarray: The calculated single span beam displacement.
            """    
            
            return (C0 * x**4)/24 - (C0 * L * x**3)/12 + ((B - A)/L + (C0 * L**3)/24) * x + A
        
        def fit_beam_displacement_two_spans(x: np.ndarray, y:np.ndarray, L0:float, p0:float=None, maxfev:int=500) -> np.ndarray:
            """ Method fits the two span beam displacement

            Arguments
            -----------
            x: np.ndarray
                distance from the selected point along the bridge
            y: np.ndarray
                displacement values
            L0: float
                length of the deck
            p0: np.ndarray, optional
                initial guess for the parameters (default is None)
            maxfev: int, optional
                maximum number of function evaluations (default is 500)
            
            Returns
            --------
            opt_params: np.ndarray
                optimized parameters for the two span beam displacement
            """
            
            opt_params, _ = curve_fit(lambda x, C0, A, B, C: two_span_beam_displacement(x, L0, C0, A, B, C), x, y, p0=p0, maxfev=maxfev)
            return opt_params
        
        def two_span_beam_displacement(x: np.ndarray, L:float, C0:float, A:float, B:float, C:float) -> np.ndarray:
            """method calculates and returns the two span beam displacement
        
            Arguments
            -----------
            x : np.ndarray
                distance from the selected point along the bridge
            L : float
                length of the deck
            C0, A, B, C : float
                interpolation variables
            
            Returns
            --------
            np.ndarray: The calculated two span beam displacement.
            """
            return np.piecewise(
                x,
                [x <= L/2, x > L/2],
                [
                    lambda x: (C0 / 24) * x**4 + (-C0 * L / 32 + 4 * (-B + (C+A)/2) / L**3) * x**3 + (C0 * L**3 / 384 + (6*B - C - 5*A) / (2*L)) * x + A,
                    lambda x: (C0 / 24) * x**4 - (C0 * (13/96) * L + 4 * (-B + (C+A)/2) / L**3) * x**3 + (C0 * (5/32) * L**2 + 12 * (-B + (C+A)/2) / L**2) * x**2 + (-C0 * (29/384) * L**3 + (9*B - 7*C/2 -11*A/2) / L) * x + C0 * (5/384) * L**4 - B + C/2 + 3*A/2
                ]
            )
        if self.data[orbit]['ndist'].size >= 3:
            if span_count == 1 :
                fit_params = fit_beam_displacement(self.data[orbit]['ndist'], self.data[orbit]['disp'], deck_length)
                solution = one_span_beam_displacement(xrange, deck_length, *fit_params)
            else:
                fit_params = fit_beam_displacement_two_spans(self.data[orbit]['ndist'], self.data[orbit]['disp'], deck_length)
                solution = two_span_beam_displacement(xrange, deck_length, *fit_params)
            
            return solution
        return None
        

    
class EW_Solver:
    """
    A class to solve the longitudinal and vertical displacement of a bridge deck using ascending and descending displacement data.
    This class initializes with time overlap information and satellite orientation and line-of-sight (LOS) information,
    and provides methods to calculate the average time series, lost longitudinal and vertical displacement,
    tilt, and deflection of the bridge deck.
    
    Attributes
    ----------
    timeOverlapInfo (dict):
        Information about the time overlap of ascending and descending data.
    theta_asc (float):
        Ascending incidence angle between LOS and nadir.
    theta_dsc (float):
        Descending incidence angle between LOS and nadir.
    alpha_asc (float):
        Ascending orbit azimuth angle.
    alpha_dsc (float):
        Descending orbit azimuth angle.
    
    Methods
    --------
    average_ts(ascending_ts: list[float], descending_ts: list[float]):
        Calculates the average time series for ascending and descending displacement data.
    ts_interpolation(combined_dates: np.ndarray, dates: np.ndarray, average_ts: np.ndarray):
        Interpolates the time series data based on combined dates and average displacement values.
    los_long_vert_displacement(interp_asc_disp: np.ndarray, interp_dsc_disp: np.ndarray, bridge_azimuth: float):
        Calculates the lost longitudinal and vertical displacement based on interpolated ascending and descending data.
    get_tilt(dataStore: dict, deck_length: float):
        Calculates the tilt of the bridge deck based on displacement data.
    get_deflection(dataStore: dict, ndist: np.ndarray, deck_length: float):
        Calculates the deflection of the bridge deck based on displacement data and normalized distances.
    """
    def __init__(
            self, 
            timeOverlapInfo : dict, 
            theta_asc : float, 
            theta_dsc : float, 
            alpha_asc : float, 
            alpha_dsc : float
            ):
        """ Initialize the EW_Solver with time overlap information and sattelite orientation and LOS information.
        
        
        Arguments
        ---------
        timeOverlapInfo : dict
            Information about the time overlap of ascending and descending data.
        theta_asc : float
            Ascending incidence angle between LOS and nadir.
        theta_dsc : float
            Descending incidence angle between LOS and nadir.
        alpha_asc : float
            Ascending orbit azimuth angle.
        alpha_dsc : float
            Descending orbit azimuth angle.
            
        """
    
        self.timeOverlapInfo = timeOverlapInfo
        self.theta_asc = np.deg2rad(theta_asc)
        self.theta_dsc = np.deg2rad(theta_dsc)
        self.alpha_asc = np.deg2rad(alpha_asc - 90)
        self.alpha_dsc = np.deg2rad(alpha_dsc - 90)
        self.combi_dates = None
        self._process_time_overlap_info()

    def _process_time_overlap_info(self) -> None:
        """ Process the time overlap information to extract ascending and descending dates and masks.
        
        This method extracts the start and end dates from the time overlap information, and creates masks for ascending and descending dates based on the specified date range. It also calculates the number of days since the start date for both ascending and descending dates.
        """

        start_date = self.timeOverlapInfo.get('rmin')
        end_date = self.timeOverlapInfo.get('rmax')
        self.__asc_time = self.timeOverlapInfo.get('ascending').get('date')
        self.__dsc_time = self.timeOverlapInfo.get('descending').get('date')
        self.__asc_mask =(self.__asc_time >= start_date) & (self.__asc_time <= end_date) 
        self.__dsc_mask = (self.__dsc_time >= start_date) & (self.__dsc_time <= end_date)
        self.__asc_num = (self.__asc_time[self.__asc_mask] - start_date).astype('timedelta64[D]').astype(int)
        self.__dsc_num = (self.__dsc_time[self.__dsc_mask] - start_date).astype('timedelta64[D]').astype(int)
        self.__combined_dates = np.sort(np.concatenate((self.__asc_num, self.__dsc_num)))
        asc_dates = self.__asc_time[self.__asc_mask]
        dsc_dates = self.__dsc_time[self.__dsc_mask]
        self.combi_dates = np.sort(np.concatenate((asc_dates, dsc_dates)))
        
    def average_ts(self, ascending_ts: list[float], descending_ts: list[float]) -> dict:
        """ Calculate the average time series for ascending and descending displacement data.
        
        This method computes the average time series based on the ascending and descending displacement data.
        
        Arguments
        ---------
        ascending_ts : list[float]
            List of ascending displacement time series data.
        descending_ts : list[float]
            List of descending displacement time series data.
        
        Returns
        --------
        dict: A dictionary containing the interpolated ascending and descending displacement time series.
        """
        ascending_ts = np.array(ascending_ts, dtype=float)
        descending_ts = np.array(descending_ts, dtype=float)
        start_date = self.timeOverlapInfo.get('rmin')
        # ascending and descending interpolated displacement
        asc_interp_disp = self.ts_interpolation(ascending_ts[self.__asc_mask], self.__asc_num,)
        dsc_interp_disp = self.ts_interpolation(descending_ts[self.__dsc_mask], self.__dsc_num,)

        if np.isnan((asc_interp_disp[0])):
            indx1 = np.where(self.__asc_time <= start_date)[0][-1]
            indx2 = np.where(self.__asc_time >= start_date)[0][0]

            x_extrap = [0, (self.__asc_time[indx2] - self.__asc_time[indx1]).days]
            y_extrap = [ascending_ts[indx1], ascending_ts[indx2]]
            extrapfunc = interp1d(x_extrap, y_extrap, kind='linear', fill_value='extrapolate')
            extrap_val = extrapfunc((start_date - self.__asc_time[indx1]).days)
            asc_interp_disp[0] = extrap_val
            #rebasing
            asc_interp_disp -= asc_interp_disp[0]

        if np.isnan((dsc_interp_disp[0])):
            indx1 = np.where(self.__dsc_time <= start_date)[0][-1]
            indx2 = np.where(self.__dsc_time >= start_date)[0][0]

            x_extrap = [0, ((self.__dsc_time[indx2] - self.__dsc_time[indx1]).days)]
            y_extrap = [descending_ts[indx1], descending_ts[indx2]]
            extrapfunc = interp1d(x_extrap, y_extrap, kind='linear', fill_value='extrapolate')
            extrap_val = extrapfunc((start_date-self.__dsc_time[indx1]).days)
            dsc_interp_disp[0] = extrap_val
            # rebasing
            dsc_interp_disp -= dsc_interp_disp[0]
        return dict(ascending = asc_interp_disp, descending = dsc_interp_disp)
            
    def ts_interpolation(self, average_ts:np.ndarray, dates: np.ndarray) -> np.ndarray:
        """ Interpolate the time series data based on combined dates and average displacement values.
        
        This method performs linear interpolation for the average time series data based on the provided dates.
        
        Arguments
        ---------
        average_ts : np.ndarray
            Average displacement time series data.
        dates : np.ndarray
            Dates corresponding to the average displacement values.
        
        Returns
        -------
        np.ndarray: Interpolated displacement values for the combined dates.
        """

        interpolated_disp = np.full_like(self.__combined_dates, np.nan, dtype=float)
        
        for i in range(len(dates) -1 ):
            # dates and displacement values of the current pair of consecutive obsevations
            x = dates[i:i+2] # consecutive dates
            y = average_ts[i:i+2] # consecutive displacement values
            # create interpolation function for the current pair
            interp_func = interp1d(x, y, kind='linear')
            #determine the indices in sorted combined_dates that fall in within this interval
            mask = (self.__combined_dates >= x[0]) & (self.__combined_dates <= x[1])
            interpolated_disp[mask] = interp_func(self.__combined_dates[mask])

        if dates[-1] < self.__combined_dates[-1]:

            x_extrap = dates[-2:] # last two dates
            y_extrap = average_ts[-2:] # last two displacement values
            extrap_func = interp1d(x_extrap, y_extrap, kind='linear', fill_value='extrapolate')
            mask_extrap = self.__combined_dates > dates[-1]
            interpolated_disp[mask_extrap] = extrap_func(self.__combined_dates[mask_extrap])

        return interpolated_disp

    def los_long_vert_displacement(self, 
                                   interp_asc_disp : np.ndarray,
                                   interp_dsc_disp : np.ndarray,  
                                   bridge_azimuth : float
                                   ) -> tuple[float, float]:
        """ 
        Calculate the lost longitudinal and vertical displacement.

        This method computes the lost displacement based on the interpolated ascending and descending data.
        
        Arguments
        ---------
        interp_asc_disp : ndarray
            Interpolated ascending displacement data.
        interp_dsc_disp : ndarray
            Interpolated descending displacement data.
        bridge_azimuth : float
            Bridge azimuth angle.
        
        Returns
        -------
        tuple: Lost longitudinal and vertical displacement.
        """
        # Convert azimuth to radians
        bridge_azimuth = np.deg2rad(bridge_azimuth)

        # Construct coefficient matrix
        A = np.array([
            [np.sin(self.theta_asc) * np.cos(self.alpha_asc - bridge_azimuth), np.cos(self.theta_asc)],
            [np.sin(self.theta_dsc) * np.cos(self.alpha_dsc - bridge_azimuth), np.cos(self.theta_dsc)]
        ])
        dLOS = np.vstack([interp_asc_disp, interp_dsc_disp])

        #solve the system of equations
        displacement = np.linalg.solve(A, dLOS)
        dL, dV = displacement[0], displacement[1]
        return dL, dV
        
    def get_tilt(self, dataStore:dict, deck_length:float) -> float:
        """ Calculate the tilt of the bridge deck.
        
        This method computes the tilt based on the vertical displacement data from the dataStore.
        
        Arguments
        ----------
        dataStore : dict
            A dictionary containing vertical displacement data for different sectors of the bridge.
        deck_length : float
            The total length of the bridge deck.
        
        Returns
        -------
            float: The tilt ratio of the bridge deck.
        """
        
        if np.isin(list(dataStore.keys()), ['N', 'S']).sum() == 2:
            tilt = abs(dataStore["N"]['vert'] - dataStore["S"]['vert']) / deck_length
            return  None if np.isnan(tilt) else tilt
        
        return None

    def get_deflection(self, dataStore:dict, ndist:np.ndarray, deck_length:float) -> float:
        """ Calculate the deflection of the bridge deck.

        This method computes the deflection based on the displacement data from the dataStore and ndist array.
        
        Arguements
        ----------
        dataStore : dict
            A dictionary containing displacement data for different sectors of the bridge.
        ndist : np.darray
            An array of normalized distances of the sectors' centroids.
        deck_length : float
            The total length of the bridge deck.
        """
        
        if np.isin(list(dataStore.keys()), ['N', "C", 'S']).sum() == 3:
            coefs = np.polyfit(ndist, [dataStore[i]['long'] for i in ['S', 'C', 'N']], deg=2)
            poly_func = np.poly1d(coefs)
            x_vals = np.linspace(ndist[0], ndist[-1], 100)
            y_vals = poly_func(x_vals)

            slope = (y_vals[-1] - y_vals[0]) / (x_vals[-1] - x_vals[0])
            intercept = y_vals[0] - slope * x_vals[0]
            deflection = y_vals - (slope * x_vals + intercept)
            max_deflec =  abs(max(deflection, key=abs)) / deck_length
            return None if np.isnan(max_deflec) else max_deflec
        return None
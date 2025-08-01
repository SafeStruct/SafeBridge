
from shapely.geometry import Polygon, LineString, Point
from numpy import ndarray
from matplotlib import pyplot 
class Plotter:
    """ 
    Plotter is a class designed to generate and customize plots for visualizing 
    bridge-related data, including persistent scatterers, deck geometry, 
    support geometries, and analytical solutions.
    Attributes:
        params (dict): A dictionary containing default plotting parameters 
            for various plot elements such as colors, line styles, markers, 
            and labels.
        _figure (matplotlib.figure.Figure): The matplotlib figure object 
            used for plotting.
        _axes (numpy.ndarray): The array of matplotlib axes objects used 
            for subplots.
    Methods:
        __init__():
            Initializes the Plotter class with default plotting parameters.
        plot(**kwargs):
            Generates plots based on the provided keyword arguments. 
            Supports different deck orientations ('NS' or 'EW') and 
            visualizes various elements such as persistent scatterers, 
            projected points, deck geometry, and analytical solutions.
            Args:
                **kwargs: Arbitrary keyword arguments containing data 
                    and configurations for the plot.
            Raises:
                ValueError: If the 'deck_orientation' argument is not 
                    'NS' or 'EW'.
        postprocess(name_tag):
            Post-processes the generated plots by setting titles, labels, 
            axis limits, and tick frequencies. Adjusts the layout and 
            appearance of the plots for better visualization.
            Args:
                name_tag (str): A string to append to plot titles for 
                    identification.
        get_figure():
            Returns the current matplotlib figure and axes objects.
            Returns:
                tuple: A tuple containing the matplotlib figure and axes 
                    objects.
    
    """
    
    def __init__(self):
        self.params = dict(
            support = dict(
                color = "darkorange", linewidth = 1.0, linestyle = "solid", alpha = 0.5, zorder = 1
            ),
            axis = dict( 
                color = "red", linewidth = 2.0, label = "Axis", linestyle = "solid", alpha = 0.2, zorder = 1
            ),
            deck = dict(
                color = "black", linewidth = 1.0, label = "Deck", linestyle = "solid", alpha = 0.5, zorder = 2,
            ),
            ascending = dict(
                color = "black", markersize = 1.5, label = "Ascending PS Points", zorder = 3,
            ),
            descending = dict(
                color = "black", markersize = 1.5, label = "Descending PS Points", zorder = 3,
            ),
            ps_graph = dict(
                color = "gray", markersize = 1.5, label = "Persistent Scatters", zorder = 3,
            ),
            projected = dict(
                color = "black", markersize = 1.5, zorder = 3,
            ),
            deck_edge = dict(
                color = "darkgreen", alpha = 0.2, zorder = 2
            ),
            sector = dict(
                color = "blue", linewidth = 1.0, alpha = 0.1, zorder = 1,
            ), 
            buffer_edge = dict(
                marker='x',  color='black',  s=50,  linewidth=2,  zorder=5,  label="Buffer intersection"
            ),
            deck_graph = dict(
                marker='+', color='black', s=70, linewidth=2, zorder=5, label="Deck intersection"
            ),
            support_graph = dict(
                marker='1', color='black', s=70, zorder=5, label='Support Intersection'
            ),
            quad_curv = dict(
                color = "blue", linewidth=2, zorder=3, label="Quadratic solution", alpha=0.9
            ),
            analytical_curv = dict(
                color='fuchsia', linewidth=2, zorder=4, label='Analytical solution', alpha=0.9
            ),
            text_info = dict(
                bbox=dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5), ha='left', va='top'
            ),
            longitudinal = dict(
              color='blue', markersize=3, label='Longitudinal Displacement', zorder=3
            ),
            vertical = dict(
              color='red', markersize=3, label='Vertical Displacement', zorder=3
            )
        )
        self._figure = None
        self._axes = None

    def plot(self, **kwargs):
        """
        Plots the data based on the provided parameters and deck orientation.
        This method generates a multi-panel plot depending on the `deck_orientation` 
        parameter. It supports two orientations: "NS" (North-South) and "EW" (East-West). 
        The plots include various geometries, projections, and analytical solutions 
        related to the deck and its components.
        Parameters:
            **kwargs: Arbitrary keyword arguments containing the data and configurations 
                      for the plot. The following keys are expected:
                - deck_orientation (str): Orientation of the deck, either "NS" or "EW".
                - descending_geom (dict): Geometry data for the descending orbit.
                - ascending_geom (dict): Geometry data for the ascending orbit.
                - descending_proj (dict): Projected geometry data for the descending orbit.
                - ascending_proj (dict): Projected geometry data for the ascending orbit.
                - axis (shapely.geometry.LineString): Longitudinal axis of the deck.
                - deck (shapely.geometry.Polygon): Geometry of the deck.
                - sectors (list): List of sector geometries.
                - support (list): List of support geometries.
                - deck_edges (list): List of deck edge geometries.
                - buf_dist (float): Buffer distance for deck edges.
                - deck_graph (list): Data for plotting the deck limit in the graph.
                - buffer_edges (list): Data for plotting buffer points in the graph.
                - support_graph (dict): Support graph data with key 'p1' for support points.
                - descending_geom_graph (dict): Graph data for descending orbit geometry.
                - ascending_geom_graph (dict): Graph data for ascending orbit geometry.
                - scaling_factor (float): Scaling factor for graph data.
                - descending_quad_solution (dict): Quadratic solution for descending orbit.
                - ascending_quad_solution (dict): Quadratic solution for ascending orbit.
                - descending_analytical_solution (list): Analytical solution for descending orbit.
                - ascending_analytical_solution (list): Analytical solution for ascending orbit.
                - descending_tilt_deflection (tuple): Tilt and deflection for descending orbit.
                - ascending_tilt_deflection (tuple): Tilt and deflection for ascending orbit.
                - timeseries (list): Time series data for longitudinal and vertical plots.
                - longitudinal (list): Longitudinal data for time series plot.
                - vertical (list): Vertical data for time series plot.
                - ew_tilt_deflection (tuple): Tilt and deflection for EW orientation.
        Raises:
            ValueError: If `deck_orientation` is not "NS" or "EW".
        Notes:
            - The method uses Matplotlib for plotting.
            - The generated figure and axes are stored in `self._figure` and `self._axes`.
            - The plot includes various elements such as scatter plots, line plots, 
              filled polygons, and text annotations.
        """
        
        
        if kwargs.get('deck_orientation', None) == "NS":
            fig, axs = pyplot.subplots(2, 2, figsize=(12, 12), dpi=300, 
                                       gridspec_kw={'height_ratios': [1, 1], 'width_ratios': [1, 1]})
        # if not generate 3 by 2 graph
        elif kwargs.get('deck_orientation', None) == "EW":
            fig, axs = pyplot.subplots(3, 2, figsize=(12, 18), 
                                       dpi=300, gridspec_kw={'height_ratios': [1, 1, 0.5], 'width_ratios': [1, 1]})
            gs = axs[0,0].get_gridspec()
            #removing the las axes from the last row
            for j in [0,1]:
                fig.delaxes(axs[2, j])
                axs[2,j] = None
            axs[2,0] = fig.add_subplot(gs[2,:])
        else:
            raise ValueError("deck_orientation must be either 'NS' or 'EW'.")
        
        # First row of the plot
        for col in range(2):
            keyword = "descending" if col == 0 else "ascending"
            # Plotting the persistent scatterers
            axs[0, col].plot(kwargs.get(keyword+"_geom")['x'],
                            kwargs.get(keyword+"_geom")['y'], 
                            "o", **self.params[keyword])
            # Plotting the projected persistent scatterers
            axs[0, col].plot(kwargs.get(keyword+"_proj")['x'],
                            kwargs.get(keyword+"_proj")['y'], 
                            "o", **self.params["projected"]) 
            # Plotting the line connecting the geom and projected points
            axs[0, col].plot([kwargs.get(keyword+"_geom")['x'],
                            kwargs.get(keyword+"_proj")['x']],
                            [kwargs.get(keyword+"_geom")['y'],
                            kwargs.get(keyword+"_proj")['y']], 
                            '--', color="black", alpha=0.2)
            # Plottting the longitudinal axis
            axs[0, col].plot(kwargs.get('axis').coords.xy[0],
                            kwargs.get('axis').coords.xy[1], **self.params["axis"])
            # Plotting the deck geometry
            axs[0, col].plot(kwargs.get('deck').exterior.xy[0],
                            kwargs.get('deck').exterior.xy[1], 
                            **self.params["deck"])
            # Plotting the sector geometries
            t = 0 
            for s in range(len(kwargs.get('sectors'))):
                axs[0, col].fill(kwargs.get('sectors')[s][0].exterior.xy[0],
                            kwargs.get('sectors')[s][0].exterior.xy[1],
                            label="Sector" if t == 0 else "",
                            **self.params["sector"])
                t += 1
            # Plotting the support geometries
            t = 0
            for sup in kwargs.get('support'):
                axs[0, col].fill(sup[0].exterior.xy[0],
                            sup[0].exterior.xy[1], 
                            label="Support" if t == 0 else "",
                            **self.params["support"])
                t += 1
            # Plotting the deck edges with buffer geometry
            for i in range(len(kwargs.get('deck_edges'))):
                axs[0, col].fill(kwargs.get('deck_edges')[i].buffer(kwargs.get('buf_dist')).exterior.xy[0],
                            kwargs.get('deck_edges')[i].buffer(kwargs.get('buf_dist')).exterior.xy[1],
                            **self.params["deck_edge"])
            
            # Second row from here

            # Plotting the deck limit in graph
            axs[1, col].scatter(kwargs.get('deck_graph'), (0,0), **self.params["deck_graph"])
            # Plotting te buffer point in graph
            axs[1, col].scatter(kwargs.get('buffer_edges'), (0,0), **self.params["buffer_edge"])
            # Plotting the supports if any
            support = kwargs.get("support_graph", {}).get('p1', ndarray(0))
            

            # plot support on the graph if it is not empty
            if support.size != 0:
                axs[1, col].scatter(support, support * 0, **self.params["support_graph"])
            orbit_key = "descending" if col == 0 else "ascending"
            # Plotting the PS graph
            axs[1, col].plot(kwargs.get(orbit_key+"_geom_graph")['x'],
                            kwargs.get(orbit_key+"_geom_graph")['y'] * kwargs.get('scaling_factor', 1),
                            "o", **self.params["ps_graph"])
            
            # Plotting the quadratic and analytical solutions if deck orientation is NS
            if kwargs.get('deck_orientation', None) == "NS":
            # Plotting the quadratic solution
                try:
                    axs[1, col].plot(kwargs.get(orbit_key+"_quad_solution")['x'][0],
                                    kwargs.get(orbit_key+"_quad_solution")['y'][0],
                                    "--", **self.params["quad_curv"])
                    axs[1, col].plot(kwargs.get(orbit_key+"_quad_solution")['x'][0],
                                    kwargs.get(orbit_key+"_analytical_solution")[0][0],
                                    '-.', **self.params["analytical_curv"])
                    # Adding tilt and deflection information
                    axs[1, col].text(0.02, 0.98,
                                f"tilt: {kwargs.get(orbit_key+'_tilt_deflection')[0]:.6f}\n"
                                f"deflection: {kwargs.get(orbit_key+'_tilt_deflection')[1]:.6f}",
                                transform=axs[1, col].transAxes, **self.params["text_info"])
                except Exception as e:
                    # logging can be added here if needed
                    pass
        
        if kwargs.get("deck_orientation", None) == "EW":
            try:
                axs[2,0].plot(kwargs.get('timeseries')[0], kwargs.get('longitudinal')[0], 
                            'o', **self.params["longitudinal"])
                axs[2,0].plot(kwargs.get('timeseries')[0], kwargs.get('vertical')[0], 
                            'x', **self.params["vertical"])
                axs[2,0].text(0.01, 0.97,
                            f"tilt: {kwargs.get('ew_tilt_deflection')[0]:.6f}\n"
                            f"deflection: {kwargs.get('ew_tilt_deflection')[1]:.6f}",
                            transform=axs[2, 0].transAxes, **self.params["text_info"])
            except Exception as e:
                # logging can be added here if needed
                pass
        self._axes = axs
        self._figure = fig
  
    def postprocess(self, name_tag):
        """Post-process the plot to set titles, labels, and limits."""
        titles = [
            f"Descending PSs along Bridge Longitudinal Axis - {name_tag}",
            f"Ascending PSs along Bridge Longitudinal Axis - {name_tag}",
            f"Descending PSs along Longitudinal Axis - {name_tag}",
            f"Ascending PSs along Longitudinal Axis - {name_tag}",
        ]
        x_labels = [
            "Longitude", 
            "Longitude",
            "Normalized distance along axis", 
            "Normalized distance along axis", 
        ]
        y_labels = [
            "Latitude", 
            "Latitude",
            "Displacement in [mm]", 
            "Displacement in [mm]", 
        ]
        
        def diff(vals):
            a,b = vals
            return abs(a-b)
        
        def mean(vals):
            a,b = vals
            return (a+b)/2
        
        # First Row
        row1_xlim = self._axes[0,0].get_xlim()
        row1_ylim = self._axes[0,0].get_ylim()

        if diff(row1_xlim) > diff(row1_ylim):
            a = mean(row1_ylim) 
            b = diff(row1_xlim) * 0.5 
            lower_limit = a - b - b * 0.1
            upper_limit = a + b + b * 0.1
            self._axes[0,0].set_ylim(lower_limit, upper_limit)
            self._axes[0,1].set_ylim(lower_limit, upper_limit)
            
            # increase the both end limits of the x and y-axes by 10 percent from the current status
            self._axes[0,0].set_xlim(row1_xlim[0] - 0.1 * diff(row1_xlim), row1_xlim[1] + 0.1 * diff(row1_xlim))
            self._axes[0,1].set_xlim(row1_xlim[0] - 0.1 * diff(row1_xlim), row1_xlim[1] + 0.1 * diff(row1_xlim))
        else:
            a = mean(row1_xlim)
            b = diff(row1_ylim) / 2
            lower_limit = a - b - b * 0.1
            upper_limit = a + b + b * 0.1
            self._axes[0,0].set_xlim(lower_limit, upper_limit)
            self._axes[0,1].set_xlim(lower_limit, upper_limit)
            self._axes[0,0].set_ylim(row1_ylim[0] - 0.1 * diff(row1_ylim), row1_ylim[1] + 0.1 * diff(row1_ylim))
            self._axes[0,1].set_ylim(row1_ylim[0] - 0.1 * diff(row1_ylim), row1_ylim[1] + 0.1 * diff(row1_ylim))
        
        # set the frequency of the ticks on the axis
        self._axes[0,0].set_xticks(self._axes[0,0].get_xticks()[::2])
        self._axes[0,1].set_xticks(self._axes[0,1].get_xticks()[::2])
              
        # Second Row
        ylims  = max([
            abs(self._axes[1,0].get_ylim()[0]),
            abs(self._axes[1,0].get_ylim()[1]),
            abs(self._axes[1,1].get_ylim()[0]),
            abs(self._axes[1,1].get_ylim()[1]),
        ])
        self._axes[1,0].set_ylim(-ylims - ylims * 0.25, ylims + ylims * 0.25)
        self._axes[1,1].set_ylim(-ylims - ylims * 0.25, ylims + ylims * 0.25)
        
        for i in range(2):
            for j in range(2):
                idx = i * 2 + j
                self._axes[i, j].set_title(titles[idx])
                self._axes[i, j].set_xlabel(x_labels[idx])
                self._axes[i, j].set_ylabel(y_labels[idx])
                self._axes[i, j].legend()

        # Third Row
        try:
            ylim = max([abs(self._axes[2,0].get_ylim()[0]),abs(self._axes[2,0].get_ylim()[1])])

            self._axes[2,0].set_ylim(-(ylim + ylim * 0.75), ylim + ylim * 0.75)
            self._axes[2,0].legend()
        
        except Exception as e:
            # this happens when there is no third row which bridge needs to be oriented in NS direction
            pass
        
    def get_figure(self):
        """Returns the current figure."""
        return (self._figure, self._axes)
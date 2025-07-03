from shapely.geometry import Polygon, LineString, Point
from numpy import ndarray
from matplotlib import pyplot 
class Plotter:
    """ A class to plot the results of the SafeBridge analysis.
    
    This class provides methods to visualize the deck geometry, sector geometries, axis geometry, support geometries, deck edges, and persistent scatterers (PS) in both ascending and descending directions. It also allows for the visualization of projected PS points, buffer edges, deck edge graphs, and support graphs. The class supports plotting of quadratic and analytical solutions for PS points, as well as tilt deflections.
    
    Attributes
    ----------
    params : dict
        A dictionary containing the parameters for plotting various elements for styling.

    _figure : matplotlib.figure.Figure
        The figure object for the plot.
    _axes : numpy.ndarray
        The axes of the plot.
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
            deck_edge_graph = dict(
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
        )
        self._figure = None
        self._axes = None

    def plot(self, 
             deck_geom,
             sector_geom,
             axis_geom,
             support_geom,
             deck_edges,
             ascending_geom,
             descending_geom,
             buf_dist,
             projected_ascending,
             projected_descending,
             buffer_edges,
             deck_edge_graph,
             ascending_geom_graph,
             descending_geom_graph,
             support_graph,
             deck_orientation,
             ascending_quad_solution,
             descending_quad_solution,
             ascending_analytic_solution,
             descending_analytic_solution,
             ascending_tilt_deflection,
             descending_tilt_deflection,
             ew_tilt_deflection
             ):
        
        fig, axs = pyplot.subplots(2, 2, figsize=(12, 12), dpi=300)

        axs[0,0].plot(descending_geom['x'],descending_geom['y'],"o", **self.params["descending"])
        axs[0,0].plot(projected_descending['x'], projected_descending['y'],"o", **self.params["projected"])
        axs[0,0].plot([descending_geom['x'], projected_descending['x']],[descending_geom['y'], projected_descending['y']], '--', color="black", alpha=0.2)
        axs[1,0].plot(descending_geom_graph['x'], descending_geom_graph['y'], "o", **self.params["ps_graph"])
        
            
        axs[0,1].plot(ascending_geom['x'],ascending_geom['y'], "o", **self.params["ascending"])
        axs[0,1].plot(projected_ascending['x'],projected_ascending['y'],"o", **self.params["projected"])
        axs[0,1].plot([ascending_geom['x'], projected_ascending['x']],[ascending_geom['y'], projected_ascending['y']], '--', color="black", alpha=0.2)
        axs[1,1].plot(ascending_geom_graph['x'], ascending_geom_graph['y'], "o", **self.params["ps_graph"])
        
        for rows in range(2):
            axs[0,rows].plot(axis_geom.coords.xy[0], axis_geom.coords.xy[1], **self.params["axis"])
            
            axs[0, rows].plot(deck_geom.exterior.xy[0], 
                              deck_geom.exterior.xy[1], 
                              **self.params["deck"],
                              )
            
            t = 0        
            for s in range(len(sector_geom)):
                axs[0, rows].fill(sector_geom[s][0].exterior.xy[0],
                            sector_geom[s][0].exterior.xy[1],
                            label="Sector" if t == 0 else "",
                            **self.params["sector"])
                t += 1

            t = 0
            for sup in support_geom:
                axs[0, rows].fill(sup[0].exterior.xy[0],
                            sup[0].exterior.xy[1], 
                            label="Support" if t == 0 else "",
                            **self.params["support"]
                            )
                t += 1
            for i in range(len(deck_edges)):
                axs[0, rows].fill(deck_edges[i].buffer(buf_dist).exterior.xy[0],
                            deck_edges[i].buffer(buf_dist).exterior.xy[1],
                            **self.params["deck_edge"])
            
            # buffer edge at the bottom
            axs[1,rows].scatter(buffer_edges, (0,0), **self.params["buffer_edge"])
            # deck edge graph at the bottom pane
            axs[1,rows].scatter(deck_edge_graph, (0,0), **self.params["deck_edge_graph"])
            # support graph geometry at the bottom pane
            if support_graph['p1'].size != 0:
                axs[1,rows].scatter(support_graph['p1'], support_graph['p1'] * 0, **self.params["support_graph"])
        
        if deck_orientation == "NS":
            try:
                axs[1,0].plot(descending_quad_solution['x'][0], descending_quad_solution['y'][0], "--",**self.params["quad_curv"])
                axs[1,0].plot(descending_quad_solution['x'][0], descending_analytic_solution, '-.', **self.params["analytical_curv"])
                axs[1,0].text(0.02, 0.98, 
                            f"tilt: {descending_tilt_deflection[0]:.6f}\ndeflection: {descending_tilt_deflection[1]:.6f}", 
                            transform=axs[1,0].transAxes, **self.params["text_info"])

                axs[1,1].plot(ascending_quad_solution['x'][0], ascending_quad_solution['y'][0], "--",**self.params["quad_curv"])
                axs[1,1].plot(ascending_quad_solution['x'][0], ascending_analytic_solution, '-.', **self.params["analytical_curv"] )
                axs[1,1].text(0.02, 0.98, 
                            f"tilt: {ascending_tilt_deflection[0]:.6f}\ndeflection: {ascending_tilt_deflection[1]:.6f}", 
                            transform=axs[1,1].transAxes, **self.params["text_info"])
            except Exception as e:
                pass
        else:
            try:
                text_tag = f"tilt: {ew_tilt_deflection[0]:.6f}\ndeflection: {ew_tilt_deflection[1]:.6f}"
                axs[1,0].text(0.02, 0.98, 
                            text_tag, 
                            transform=axs[1,0].transAxes, **self.params["text_info"])
            except Exception as e:
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
        ylim = max([
            abs(self._axes[1,0].get_ylim()[0]),
            abs(self._axes[1,0].get_ylim()[1]),
            abs(self._axes[1,1].get_ylim()[0]),
            abs(self._axes[1,1].get_ylim()[1]),
        ])

        for i in range(2):
            for j in range(2):
                idx = i * 2 + j
                self._axes[i, j].set_title(titles[idx])
                self._axes[i, j].set_xlabel(x_labels[idx])
                self._axes[i, j].set_ylabel(y_labels[idx])
                self._axes[i, j].legend()

                if (i,j) in [(0,0), (0,1)]:
                    self._axes[i,j].set_xlim(self._axes[i,j].get_xlim()[0] - 10 , self._axes[i,j].get_xlim()[1] + 10)
                    self._axes[i,j].set_ylim(self._axes[i,j].get_ylim()[0] - 10 , self._axes[i,j].get_ylim()[1] + 10)
                    self._axes[i,j].ticklabel_format(useOffset=False, style='plain')
                    

                    self._axes[i, j].axes.set_aspect('equal')

                if (i,j) in [(1,0), (1,1)]:
                    self._axes[i, j].axes.set_ylim([-ylim-5, ylim+5])
        
        pyplot.tight_layout()
        
    def get_figure(self):
        """Returns the current figure."""
        return (self._figure, self._axes)
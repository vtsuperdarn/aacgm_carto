import cartopy
from cartopy.mpl.geoaxes import GeoAxes
import aacgmv2
import numpy
from shapely.geometry import  MultiLineString, mapping, LineString, Polygon
from matplotlib.projections import register_projection
import copy
import datetime


class AxesAACGM(GeoAxes):
    name = 'aacgmv2'
        
    def __init__(self, *args, **kwargs):
        if "map_projection" in kwargs:
            map_projection = kwargs.pop("map_projection")
        else:
            map_projection = cartopy.crs.NorthPolarStereo()
            print("map_projection keyword not set, setting it to cartopy.crs.NorthPolarStereo()")
        # first check if datetime keyword is given!
        # it should be since we need it for aacgm
        if "plot_date" in kwargs:
            self.plot_date = kwargs.pop("plot_date")
        else:
            raise TypeError("need to provide a date using 'plot_date' keyword for aacgmv2 plotting")
        # Now work with the coords!
        supported_coords = [ "geo", "aacgmv2", "aacgmv2_mlt" ]
        if "coords" in kwargs:
            self.coords = kwargs.pop("coords")
            if self.coords not in supported_coords:
                err_str = "coordinates not supported, choose from : "
                for _n,_sc in enumerate(supported_coords):
                    if _n + 1 != len(supported_coords):
                        err_str += _sc + ", "
                    else:
                        err_str += _sc
                raise TypeError(err_str)
        else:
            self.coords = "aacgmv2"
            print("coords keyword not set, setting it to aacgmv2")
        # finally, initialize te GeoAxes object
        GeoAxes.__init__(self, map_projection=map_projection,*args, **kwargs)

    def overaly_coast_lakes(self, resolution='110m', color='black', **kwargs):
        """
        Overlay AACGM coastlines and lakes
        """
        kwargs['edgecolor'] = color
        kwargs['facecolor'] = 'none'
        # overaly coastlines
        feature = cartopy.feature.NaturalEarthFeature('physical', 'coastline',
                                                      resolution, **kwargs)
        self.add_feature( cartopy.feature.COASTLINE, **kwargs )
        self.add_feature( cartopy.feature.LAKES, **kwargs )
        
    def coastlines(self,resolution='110m', color='black', **kwargs):
        # details!
        kwargs['edgecolor'] = color
        kwargs['facecolor'] = 'none'
        feature = cartopy.feature.NaturalEarthFeature('physical', 'coastline',
                                                      resolution, **kwargs)
        return self.add_feature(feature, **kwargs)
        
    def add_feature(self, feature, **kwargs):
        aacgm_geom = self.get_aacgm_geom(feature)
        aacgm_feature = cartopy.feature.ShapelyFeature(aacgm_geom, cartopy.crs.Geodetic(),\
                                                 **kwargs)
        # Now we'll set facecolor as None because aacgm doesn't close
        # continents near equator and it turns into a problem
        if 'edgecolor' not in kwargs:
            kwargs['edgecolor'] = "black"
        if 'facecolor' in kwargs:
            print("manually setting facecolor keyword to none as aacgm fails for fill! want to know why?? think about equator!")
        kwargs['facecolor'] = "none"
        super().add_feature(aacgm_feature, **kwargs)
        
    def get_aacgm_geom(self, feature, out_height=300. ):
        new_i = []
        # cartopy.feature.COASTLINE
        for _n,i in enumerate(feature.geometries()):
            aa = mapping(i)
            mag_list = []
            geo_coords = aa["coordinates"]
            for _ni, _list in enumerate(geo_coords):
                mlon_check_jump_list = []
                split_mag_list = None
                if len(_list) == 1:
                    _loop_list = _list[0]
                else:
                    _loop_list = _list
                for _ngc, _gc in enumerate(_loop_list):
                    _mc = aacgmv2.get_aacgm_coord(_gc[1], _gc[0], out_height, self.plot_date)
                    if numpy.isnan(_mc[0]):
                        continue 
                    mlon_check_jump_list.append( _mc[1] )

                    mag_list.append( (_mc[1], _mc[0]) )
                # check for unwanted jumps
                mlon_check_jump_list = numpy.array( mlon_check_jump_list )

                jump_arr = numpy.diff( mlon_check_jump_list )
                bad_inds = numpy.where( numpy.abs(jump_arr) > 10.)[0]
                # delete the range of bad values
                # This is further complicated because
                # in some locations mlon jumps from -177 to +178
                # and this causes jumps in the maps! To deal with 
                # this we'll split arrays of such jumps 
                # (these jumps typically have just one bad ind )
                # and make them into two seperate entities (LineStrings)
                # so that shapely will treat them as two seperate boundaries!
                if len(bad_inds) > 0:
                    if len(bad_inds) > 1:
                        mag_list = [i for j, i in enumerate(mag_list) if j-1 not in numpy.arange(bad_inds[0], bad_inds[1])]
                    else:
                        split_mag_list = mag_list[bad_inds[0]+1:]
                        mag_list = mag_list[:bad_inds[0]+1]
                mag_coords = tuple(mag_list)
                if len(mag_list) > 1:
                    new_i.append( mag_coords )
                if split_mag_list is not None:
        #             print(split_mag_list)
                    if len(split_mag_list) > 1:
                        new_i.append( tuple(split_mag_list) )
        aacgm_coast = MultiLineString( new_i )
        return aacgm_coast

    def mark_latitudes(self, lat_arr, lon_location=45, **kwargs):
        """
        mark the latitudes
        Write down the latitudes on the map for labeling!
        we are using this because cartopy doesn't have a 
        label by default for non-rectangular projections!
        """
        if isinstance(lat_arr, list):
            lat_arr = numpy.array(lat_arr)
        else:
            if not isinstance(lat_arr, numpy.ndarray):
                raise TypeError('lat_arr must either be a list or numpy array')
        # make an array of lon_location
        lon_location_arr = numpy.full( lat_arr.shape, lon_location )
        proj_xyz = self.projection.transform_points(\
                            cartopy.crs.Geodetic(),\
                            lon_location_arr,\
                            lat_arr
                            )
        # plot the lats now!
        out_extent_lats = False
        for _np,_pro in enumerate(proj_xyz[..., :2].tolist()):
            # check if lats are out of extent! if so ignore them
            lat_lim = self.get_extent(crs=cartopy.crs.Geodetic())[2::]
            if (lat_arr[_np] >= min(lat_lim)) and (lat_arr[_np] <= max(lat_lim)):
                self.text( _pro[0], _pro[1], str(lat_arr[_np]), **kwargs)
            else:
                out_extent_lats = True
        if out_extent_lats:
            print( "some lats were out of extent ignored them" )

    def mark_longitudes(self, lon_arr=numpy.arange(-180,180,60), **kwargs):
        """
        mark the longitudes
        Write down the longitudes on the map for labeling!
        we are using this because cartopy doesn't have a 
        label by default for non-rectangular projections!
        This is also trickier compared to latitudes!
        """
        if isinstance(lon_arr, list):
            lon_arr = numpy.array(lon_arr)
        else:
            if not isinstance(lon_arr, numpy.ndarray):
                raise TypeError('lat_arr must either be a list or numpy array')

        # get the boundaries
        [x1, y1], [x2, y2] = self.viewLim.get_points()
        bound_lim_arr = []
        right_bound = LineString(([-x1, y1], [x2, y2]))
        top_bound = LineString(([x1, -y1], [x2, y2]))
        bottom_bound = LineString(([x1, y1], [x2, -y2]))
        left_bound = LineString(([x1, y1], [-x2, y2]))
        plot_outline = MultiLineString( [\
                                        right_bound,\
                                        top_bound,\
                                        bottom_bound,\
                                        left_bound\
                                        ] )

        # bound_lim_arr.append( LineString(([-x1, y1], [x2, y2])) )# right
        # bound_lim_arr.append( LineString(([x1, -y1], [x2, y2])) )#top
        # bound_lim_arr.append( LineString(([x1, y1], [x2, -y2])) )# bottom
        # bound_lim_arr.append( LineString(([x1, y1], [-x2, y2])) )# left
        # plot_outline = MultiLineString( bound_lim_arr )
        # get the plot extent, we'll get an intersection
        # to locate the ticks!
        plot_extent = self.get_extent(cartopy.crs.Geodetic())
        line_constructor = lambda t, n, b: numpy.vstack(\
                        (numpy.zeros(n) + t, numpy.linspace(b[2], b[3], n))\
                        ).T
        for t in lon_arr:
            xy = line_constructor(t, 30, plot_extent)
            # print(xy)
            proj_xyz = self.projection.transform_points(\
                            cartopy.crs.PlateCarree(), xy[:, 0], xy[:, 1]\
                            )
            xyt = proj_xyz[..., :2]
            ls = LineString(xyt.tolist())
            locs = plot_outline.intersection(ls)
            if not locs:
                continue
            # we need to get the alignment right
            # so get the boundary closest to the label
            # and plot it!
            closest_bound =min( [\
                            right_bound.distance(locs),\
                            top_bound.distance(locs),\
                            bottom_bound.distance(locs),\
                            left_bound.distance(locs)\
                            ] )
            if closest_bound == right_bound.distance(locs):
                ha = 'left'
                va = 'top'
            elif closest_bound == top_bound.distance(locs):
                ha = 'left'
                va = 'bottom'
            elif closest_bound == bottom_bound.distance(locs):
                ha = 'left'
                va = 'top'
            else:
                ha = 'right'
                va = 'top'
            self.text( locs.bounds[0],locs.bounds[1], str(t), ha=ha, va=va)


        
# Now register the projection with matplotlib so the user can select
# it.
register_projection(AxesAACGM)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import datetime
    import cartopy
    import numpy

    plot_date = datetime.datetime(2015,1,1)
    fig = plt.figure()
    ax = fig.add_subplot(
                projection='aacgmv2',\
                map_projection = cartopy.crs.NorthPolarStereo(),\
                coords="aacgmv2", plot_date=plot_date
                )
    # uncomment lines below to add coastlines and lakes individually
    # ax.coastlines()
    # ax.add_feature( cartopy.feature.LAKES)
    # or add coastlines and lakes together!
    ax.overaly_coast_lakes()
    # plot set the map bounds
    ax.set_extent([-180, 180, 40, 90], crs=cartopy.crs.PlateCarree())
    # plot a random line!
    # ax.scatter(54, 60, transform=cartopy.crs.Geodetic())
    ax.plot( [-175, 175], [60,60], transform=cartopy.crs.Geodetic() )
    # overaly gridlines!
    ax.gridlines(linewidth=0.5)
    ax.mark_latitudes(numpy.arange(20,90,10), fontsize=10)
    ax.mark_longitudes(fontsize=10)
    plt.savefig("test_plots/carto_test.png")
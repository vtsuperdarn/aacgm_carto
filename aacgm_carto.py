import cartopy
from cartopy.mpl.geoaxes import GeoAxes
import aacgmv2
import numpy
from shapely.geometry import  MultiLineString, mapping, shape
from matplotlib.projections import register_projection
import copy
import datetime


class AxesAACGM(GeoAxes):
    name = 'aacgmv2'
        
    def __init__(self, *args, **kwargs):
        if "map_projection" in kwargs:
            map_projection=kwargs.pop("map_projection")
        else:
            map_projection = cartopy.crs.NorthPolarStereo()
            print("map_projection keyword not set, setting it to cartopy.crs.NorthPolarStereo()")
        GeoAxes.__init__(self, map_projection=map_projection,*args, **kwargs)

    def overaly_coast_lakes(self, coords="aacgmv2", plot_date=None,\
                  resolution='110m', color='black', **kwargs):
        """
        Overlay AACGM coastlines and lakes
        """
        supported_coords = [ "geo", "aacgmv2", "aacgmv2_mlt" ]
        if coords not in supported_coords:
            print("coordinates not supported, choose from ", supported_coords)
            return
        if (coords != "geo") and (datetime == None):
            print("Your need to set datetime keyword for plotting in AACGM")
            return

        kwargs['edgecolor'] = color
        kwargs['facecolor'] = 'none'
        # overaly coastlines
        feature = cartopy.feature.NaturalEarthFeature('physical', 'coastline',
                                                      resolution, **kwargs)
        self.add_feature( cartopy.feature.COASTLINE, coords=coords,\
                                plot_date=plot_date, **kwargs )
        self.add_feature( cartopy.feature.LAKES, coords=coords,\
                                plot_date=plot_date, **kwargs )
        
    def coastlines(self, coords="aacgmv2", plot_date=None,\
                  resolution='110m', color='black', **kwargs):
        supported_coords = [ "geo", "aacgmv2", "aacgmv2_mlt" ]
        if coords not in supported_coords:
            print("coordinates not supported, choose from ", supported_coords)
            return
        if (coords != "geo") and (datetime == None):
            print("Your need to set datetime keyword for plotting in AACGM")
            return
        # actual details!
        kwargs['edgecolor'] = color
        kwargs['facecolor'] = 'none'
        feature = cartopy.feature.NaturalEarthFeature('physical', 'coastline',
                                                      resolution, **kwargs)
        return self.add_feature(feature, coords=coords,\
                                plot_date=plot_date, **kwargs)
        
    def add_feature(self, feature, coords="aacgmv2", plot_date=None, **kwargs):
        supported_coords = [ "geo", "aacgmv2", "aacgmv2_mlt" ]
        if coords not in supported_coords:
            print("coordinates not supported, choose from ", supported_coords)
            return
        if (coords != "geo") and (datetime == None):
            print("Your need to set datetime keyword for plotting in AACGM")
            return
        aacgm_geom = self.get_aacgm_geom(feature, plot_date)
        aacgm_feature = cartopy.feature.ShapelyFeature(aacgm_geom, cartopy.crs.PlateCarree(),\
                                                 **kwargs)
        # Now we'll set facecolor as None because aacgm doesn't close
        # continents near equator and it turns into a problem
        if 'edgecolor' not in kwargs:
            kwargs['edgecolor'] = "black"
        if 'facecolor' in kwargs:
            print("manually setting facecolor keyword to none as aacgm fails for fill! want to know why?? think about equator!")
        kwargs['facecolor'] = "none"
        super().add_feature(aacgm_feature, **kwargs)
        
    def get_aacgm_geom(self, feature, plot_date, out_height=300. ):
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
                    _mc = aacgmv2.get_aacgm_coord(_gc[1], _gc[0], out_height, plot_date)
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

        
# Now register the projection with matplotlib so the user can select
# it.
register_projection(AxesAACGM)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import datetime
    import cartopy

    plot_date = datetime.datetime(2015,1,1)
    fig = plt.figure()
    ax = fig.add_subplot(projection='aacgmv2',map_projection = cartopy.crs.NorthPolarStereo())
    # uncomment lines below to add coastlines and lakes individually
    # ax.coastlines(coords="aacgmv2", plot_date=plot_date)
    # ax.add_feature( cartopy.feature.LAKES, plot_date=plot_date )
    # or add coastlines and lakes together!
    ax.overaly_coast_lakes(coords="aacgmv2", plot_date=plot_date)
    # plot set the map bounds
    ax.set_extent([-180, 180, 40, 90], crs=cartopy.crs.PlateCarree())
    # plot a random line!
    ax.plot([-0.08, 132], [51.53, 43.17], transform=cartopy.crs.PlateCarree())
    # overaly gridlines!
    ax.gridlines(linewidth=0.5)
    plt.savefig( "test_plots/carto_test.png" )
VECTORIZED HAVERSINE METHOD TO CALCULATE SPACE ANGLE GIVEN (RA, DEC) OF TWO POINTS
def hvvec(lon1, lat1, lon2, lat2):

    #Convert decimal degrees to Radians:
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)

    #Implementing Haversine Formula: 
    dlon = np.subtract(lon2, lon1)
    dlat = np.subtract(lat2, lat1)

    a = np.add(np.power(np.sin(np.divide(dlat, 2)), 2),  
                          np.multiply(np.cos(lat1), 
                                      np.multiply(np.cos(lat2), 
                                                  np.power(np.sin(np.divide(dlon, 2)), 2))))
    c = np.multiply(2, np.arcsin(np.sqrt(a)))

    return c



ASTROPY METHOD

asp = astropy.coordinate.position_angle( np.deg2rad(msra), np.deg2rad(msdec), np.deg2rad(icra), np.deg2rad(icdec))

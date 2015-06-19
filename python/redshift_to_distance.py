def redshift_to_distance_integrand ( y, OmegaM, OmegaL ) :
    return ( 1 / ( y*(OmegaM/y+OmegaL*y**2)**0.5 ) )

def redshift_to_distance (Z,  h = 1, OmegaM = 0.25, OmegaL = 0.75) :
    from scipy.integrate import quad
    z = np.asarray( Z )
    if z.size == 1 : z = np.asarray( [Z] )
    r = np.zeros( z.size )
    for i in range( 0, z.size ) :
        r[i] = quad( redshift_to_distance_integrand, 1./(1.+z[i]), 1., args=(OmegaM, OmegaL) )[0] * 3000 / h
    return r

def visible_to_absolute ( m_, r_, z_, k_=0 ) :
    m = np.asarray( m_ )
    r = np.asarray( r_ )
    z = np.asarray( z_ )
    k = np.asarray( k_ )
    if ( m.size != r.size ) | ( m.size != z.size ) | ( ( k.size > 1) & ( k.size != m.size ) ) :
        print( "Wrong size of arguments" )
        raise RuntimeError
    if ( m.size == 1 ) :
        m = np.asarray( [m_] )
        r = np.asarray( [r_] )
        z = np.asarray( [z_] )
        k = np.asarray( [k_] )
    M = np.empty( m.size )
    for i in range( 0, m.size ) :
        if r[i] <= 0 :
            M[i] = 0 
        else :
            M[i] = m[i] - 5 * np.log10( r[i] * (1 + z[i]) ) - 25 
    if k.any != 0 :
        M = M - k
    return M 


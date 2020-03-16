'''
here the properties from the vterm paper are stored and can be returned as a dict
'''

from IPython.core.debugger import Tracer; debug = Tracer()


#define particle class
def init_class():
    class particle(object):
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)
#define objects for particle types
    p = dict()
    p["column"] = particle(      nu_SB	=  0.,
            mu_SB       = 0.333,
            a_geo       =  (1./0.046)**(1./2.07),
            b_geo       =  1./2.07,
            xmax        =  1.00e-05, #& !..x_max..maximale Teilchenmasse D=???e-2m
            xmin        =  1.00e-12, #& !..x_min..minimale Teilchenmasse D=200e-6m
            a_Atlas     =1.629,
            b_Atlas     =1.667,
            c_Atlas     =1586.0)
    p["plate"] = particle(      nu_SB	=  0.,
            mu_SB       = 0.3333,
            a_geo       =  (1./0.788)**(1./2.48),
            b_geo       =  1./2.48,
            xmax        =  1.00e-05, #& !..x_max..maximale Teilchenmasse D=???e-2m
            xmin        =  1.00e-12, #& !..x_min..minimale Teilchenmasse D=200e-6m
            a_Atlas     =2.265,
            b_Atlas     =2.275,
            c_Atlas     =771.138)
    p["dendrite"] = particle(    nu_SB	=  0.,
            mu_SB       = 0.3333,
            a_geo       =  (1./0.074)**(1./2.33),
            b_geo       =  1./2.33,
            xmax        =  1.00e-05, #& !..x_max..maximale Teilchenmasse D=???e-2m
            xmin        =  1.00e-12, #& !..x_min..minimale Teilchenmasse D=200e-6m
            a_Atlas     =1.133,
            b_Atlas     =1.153,
            c_Atlas     =1177.0)
    p["needle"] = particle(      nu_SB       =  0.,
            mu_SB       = 0.3333,
            a_geo       =  (1./0.005)**(1./1.89),
            b_geo       =  1./1.89,
            xmax        =  1.00e-05, #& !..x_max..maximale Teilchenmasse D=???e-2m
            xmin        =  1.00e-12, #& !..x_min..minimale Teilchenmasse D=200e-6m
            a_Atlas     =0.848,
            b_Atlas     =0.871,
            c_Atlas     =2278.0)
    '''
    = particle(       nu_SB	=  0.,
            mu_SB       = 0.3333,
            a_geo       =  (1./)**(1./),
            b_geo       =  1./,
            xmax        =  1.00e-05, #& !..x_max..maximale Teilchenmasse D=???e-2m
            xmin        =  1.00e-12, #& !..x_min..minimale Teilchenmasse D=200e-6m
            a_Atlas     =,
            b_Atlas     =,
            c_Atlas     =)
    '''
#,
    p["aggPlate"] = particle(    nu_SB     =  0.,
            mu_SB       = 0.5,
            a_geo       =  (1./0.076)**(1./2.22),
            b_geo       =  1./2.22,
            xmax        =  2.00e-05, #& !..x_max..minimale Teilchenmasse 
            xmin        =  1.00e-10, #& !..x_min..minimale Teilchenmasse
            a_Atlas     =1.366,
            b_Atlas     =1.391,
            c_Atlas     =1285.6,
            ssrg        ="0.18_0.89_2.06_0.08")
    p["aggDendrite"] = particle( nu_SB     =  0.,
            mu_SB       = 0.5,
            a_geo       =  (1./0.027)**(1./2.22),
            b_geo       =  1./2.22,
            xmax        =  2.00e-05, #& !..x_max..minimale Teilchenmasse 
            xmin        =  1.00e-10, #& !..x_min..minimale Teilchenmasse
            a_Atlas     =0.880,
            b_Atlas     =0.895,
            c_Atlas     =1393.0,
            ssrg        ="0.23_0.75_1.88_0.10")
    p["aggNeedle"] = particle(         nu_SB     =  0.,
            mu_SB       = 0.5,
            a_geo       =  (1./0.028)**(1./2.11),
            b_geo       =  1./2.11,
            xmax        =  2.00e-05, #& !..x_max..minimale Teilchenmasse 
            xmin        =  1.00e-10, #& !..x_min..minimale Teilchenmasse
            a_Atlas     =1.118,
            b_Atlas     =1.133,
            c_Atlas     =1659.5,
            ssrg        ="0.25_0.76_1.66_0.04")
    p["aggColumn"] = particle(   nu_SB     =  0.,
            mu_SB       = 0.5,
            a_geo       =  (1./0.074)**(1./2.15),
            b_geo       =  1./2.15,
            xmax        =  2.00e-05, #& !..x_max..minimale Teilchenmasse 
            xmin        =  1.00e-10, #& !..x_min..minimale Teilchenmasse
            a_Atlas     =1.583,
            b_Atlas     =1.600,
            c_Atlas     =1419.2,
            ssrg        ="0.23_1.45_2.05_0.02")
    p["Mix1"] = particle(     nu_SB     =  0.,
            mu_SB       = 0.5,
            a_geo       =  (1./0.045)**(1./2.16),
            b_geo       =  1./2.16,
            xmax        =  2.00e-05, #& !..x_max..minimale Teilchenmasse 
            xmin        =  1.00e-10, #& !..x_min..minimale Teilchenmasse
            a_Atlas     =1.233,
            b_Atlas     =1.250,
            c_Atlas     =1509.5)
    p["Mix2"] = particle(      nu_SB	=  0.,
            mu_SB       = 0.5,
            a_geo       =  (1./0.017)**(1./1.95),
            b_geo       =  1./1.95,
            a_A         =  0.066     , #projected area A=a_A*D**b_A
            b_A         =  1.79      , #projected area A=a_A*D**b_A
            xmax        =  2.00e-05, #& !..x_max..minimale Teilchenmasse 
            xmin        =  1.00e-10, #& !..x_min..minimale Teilchenmasse
            a_Atlas     =1.121,
            b_Atlas     =1.119,
            c_Atlas     =2292.2,
            ssrg        ="0.22_0.60_1.81_0.11")
    p["Mix2SB"] = particle(      nu_SB	=  0., #Mix2 ssrg and Atlas fit but SB m-D
            mu_SB       = 0.5,
            a_geo       =  5.13, #(1./0.017)**(1./1.95),
            b_geo       =  0.5,
            xmax        =  2.00e-05, #& !..x_max..minimale Teilchenmasse 
            xmin        =  1.00e-10, #& !..x_min..minimale Teilchenmasse
            a_Atlas     =1.121,
            b_Atlas     =1.119,
            c_Atlas     =2292.2,
            ssrg        ="0.23_0.60_1.8_0.11")
    p["snowSBB"] = particle(      nu_SB	=  0., #Mix2 ssrg and Atlas fit but SB m-D
            mu_SB       = 0.5,
            a_geo       =  5.13, #(1./0.017)**(1./1.95),
            b_geo       =  0.5,
            xmax        =  2.00e-05, #& !..x_max..minimale Teilchenmasse 
            xmin        =  1.00e-10, #& !..x_min..minimale Teilchenmasse
            a_vel       =   8.294,
            b_vel       =   0.125)
    '''
    p["agg"] = particle(         nu_SB     =  0.,
            mu_SB       = 0.333,
            a_geo       =  (1./)**(1./),
            b_geo       =  1./,
            xmax        =  1.00e-05, #& !..x_max..maximale Teilchenmasse D=???e-2m
            xmin        =  1.00e-12, #& !..x_min..minimale Teilchenmasse D=200e-6m
            a_Atlas     =,
            b_Atlas     =,
            c_Atlas     =)
    '''
#convert from N(m) to N(D) space
    for key in p.keys():
        curr_cat=p[key]
#for curr_cat in [cloud_ice,snow]:
        curr_cat.a_ms = (1./curr_cat.a_geo)**(1./curr_cat.b_geo)
        curr_cat.b_ms = 1./curr_cat.b_geo
        curr_cat.mu =  curr_cat.b_ms*curr_cat.nu_SB+curr_cat.b_ms-1
        curr_cat.gam = curr_cat.b_ms*curr_cat.mu_SB
    return p

def conv_from_Nm_to_ND(curr_cat,nu_SB=None):
   '''
   convert from N(m) to N(D)
   ARGUMENTS:
   curr_cat: a "particle"-object (e.g. p["Mix2"])
   nu_SB: modification of the nu_SB parameter
   '''
   if nu_SB!=None:
       curr_cat.nu_SB=nu_SB
   #convert from N(m) to N(D) space
   curr_cat.a_ms = (1./curr_cat.a_geo)**(1./curr_cat.b_geo)
   curr_cat.b_ms = 1./curr_cat.b_geo
   curr_cat.mu =  curr_cat.b_ms*curr_cat.nu_SB+curr_cat.b_ms-1
   curr_cat.gam = curr_cat.b_ms*curr_cat.mu_SB
   return curr_cat



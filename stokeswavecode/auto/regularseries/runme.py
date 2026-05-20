import auto as a
import numpy as np

# import code
#     code.interact(local=locals())

f = 0.3


sysname  = "RegularStokes"
# unames   = dict( enumerate( ["U"], start=1) )
parnames = dict( enumerate( ["ep", "mu", "B"], start=1) )

prob = a.load(
    e = sysname,
    # unames = unames,
    parnames = parnames,
    #dat='wrinkle.dat',
    #dat='start_tau20_n8_ep0.17036.dat',
    NDIM=  81, IPS = 0, IRS = 0, ILP = 0,
    ICP =  ["ep", "mu", "B"],
    NTST=  50, NCOL= 4, IAD = 3, ISP = 2, ISW = 1, IPLT= 0, NBC= 5, NINT= 1,
    NMX=  100, RL0 = 1.e-8, RL1 = 1.e+0, A0 = -1.e+8, A1 = 1.e+8, 
    NPR=  10, MXBF = 10, IID = 4, ITMX= 8, ITNW= 5, NWTN= 3, JAC= 0,
    EPSL= 1.e-6, EPSU = 1.e-6, EPSS = 1.e-5,
    DS  =  0.05, DSMIN=  0.001, DSMAX= 0.1, IADS= 1,
    NPAR = 3, THL = {}, THU = {},
    )

# Vary nu
rtot = a.run(prob)
# r2 = a.run(r1(2), DS=-0.01, DSMIN=0.001, DSMAX=0.01)
# # Use this for backwards on the 0.38 special
# # r2 = a.run(prob, DS=-0.01, DSMIN=0.001, DSMAX=0.01)
# rtot1 = a.merge(r1+r2)

# # Vary E
# r3 = a.run(prob, ICP=["E", "nu"], RL1=1.e+2, NMX=100)
# r4 = a.run(prob, ICP=["E", "nu"], RL1=1.e+2, DS=-0.1, NMX=1000)
# rtot2 = a.merge(r3+r4)
# rtot = a.merge(rtot1 + rtot2)

# nu1 = rtot1["nu"]
# E1 = rtot1["E"]
# S11 = rtot1["S1"]
# S21 = rtot1["S2"]
# S31 = rtot1["S3"]

# myfile = "auto_f"+"{:.3f}".format(f) + "_varynu.dat"
# np.savetxt(myfile, np.c_[nu1, E1, S11, S21, S31], header='[nu, E, S1, S2, S3]', comments='#')

# nu2 = rtot2["nu"]
# E2 = rtot2["E"]
# S12 = rtot2["S1"]
# S22 = rtot2["S2"]
# S32 = rtot2["S3"]
# myfile = "auto_f"+"{:.3f}".format(f) + "_varyE.dat"
# np.savetxt(myfile, np.c_[nu2, E2, S12, S22, S32], header='[nu, E, S1, S2, S3]', comments='#')

# from pylab import *
# ion()
# plot(rtot["nu"], rtot["S1"])
# plot(rtot["nu"], rtot["S2"])
# plot(rtot["nu"], rtot["S3"])

# Output columns
# savetxt('auto_f0.3_3d_A.dat', c_[nu, E, Ssq, S1, S2, S3])
# savetxt('auto_f0.3_3d_A.dat', c_[nu, E, Ssq, S1, S2, S3])
# a.plot(rtot, use_labels=False, bifurcation_x='nu', bifurcation_y='E', bifurcation_z='S2')
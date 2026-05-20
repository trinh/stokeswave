import auto as a

sysname  = "wrinkle"
unames   = dict( enumerate( ["W", "Wp", "Wpp", "Wppp", "x"], start=1) )
parnames = dict( enumerate( ["eps", "n", "tau", "Bo", "L", "p"], start=1) )

prob = a.load(
    e = sysname,
    unames = unames,
    parnames = parnames,
    dat='wrinkle.dat',
    NDIM=   len(unames), IPS = 4, IRS = 0, ILP = 1,
    ICP =  ["eps", "n", "tau"],
    NTST=  100, NCOL= 4, IAD = 3, ISP = 3, ISW = 1, IPLT= 0, NBC= 5, NINT= 1,
    NMX=  100, RL0 = 0, RL1 = 1.e+8, A0 = -1.e+8, A1 = 1.e+8, 
    NPR=  10, MXBF = 10, IID = 2, ITMX= 8, ITNW= 5, NWTN= 3, JAC= 0,
    EPSL= 1.e-8, EPSU = 1.e-8, EPSS = 1.e-6,
    DS  =  0.01, DSMIN=  0.001, DSMAX= 0.1, IADS= 1,
    NPAR = 6, THL = {}, THU = {},
    )

# # We first run the code backwards to detect the first fold point
# r1 = a.run(prob, SP=['LP1'])

# # Beginning from fold point, solve downwards in tau
# rtau = a.run(r1, IRS='LP1', ISW=-1, SP=['LP'], ICP=["tau", "eps", "n"], RL0=8, DS=-0.1, DSMAX=0.2, LAB=1)

rstart = a.run(prob, ICP=["eps", "n", "tau"])
rstart = rstart + a.run(prob, ICP=["eps", "n", "tau"], DS=-0.01)
rtot = rstart

rlast = rstart
for j in range(1,100):

    r2 = a.run(rlast, SP=['BP0'], IRS='LP1', ICP=["eps", "tau", "n"], ISW=2, DS=-0.01, DSMAX=0.5)
    rtot = rtot + r2

    # r3 = a.run(r2, SP='LP1', ICP=["eps", "n", "tau"], ISW=2, DS=0.01)
    # rtot = rtot + r3
    
    r1 = a.run(r2, SP=['BP0'], ICP=["n", "eps", "tau"], NMX = 50, ISW=-1, DS=0.01, DSMAX=0.5)
    rb = a.run(r2, SP=['BP0'], ICP=["n", "eps", "tau"], NMX = 50, ISW=-1, DS=-0.01, DSMAX=0.5)
    r1 = a.run(r1, SP = ['BP0', 'LP1'], ICP=["eps", "n", "tau"], NMX = 1000, ISW=-1, DS=0.1, DSMAX=0.1)
    # r1 = r1 + a.run(r2, ICP=["eps", "n", "tau"], ISW=-1, DS=-0.01)  
    rtot = rtot + r1 + rb
    rlast = r1

    if j == 1:
        nmat = [r1('LP1')['n']]
        epsmat = [r1('LP1')['eps']]
        taumat = [r1('LP1')['tau']]
    else:
        nmat.append(r1('LP1')['n'])
        epsmat.append(r1('LP1')['eps'])
        taumat.append(r1('LP1')['tau'])

    # r2 = a.run(r1, IRS='LP1', ICP=["eps", "tau", "n"], ISW=2, DS=-0.1)
    # rtot = rtot + r2

    # r1 = a.run(r2, ICP=["n", "eps", "tau"], NMX = 20, ISW=-1, DS=0.01)
    # r1 = a.run(r1, SP = ['LP1'], ICP=["eps", "n", "tau"], NMX = 100, ISW=-1, DS=0.1)
    # # r1 = r1 + a.run(r2, ICP=["eps", "n", "tau"], ISW=-1, DS=-0.01)  
    # rtot = rtot + r1

    # r2 = a.run(r1, IRS='LP1', ICP=["eps", "tau", "n"], ISW=2, DS=-0.1)
    # rtot = rtot + r2

    # r4 = a.run(r2, SP=1, ICP=["eps", "n", "tau"], ISW=2, DS=-0.01)
    # # r2 = r2 + a.run(solution, IRS='LP1', ICP=["eps", "tau", "n"], ISW=2, DS=-0.01)
    # rtot = rtot + r4

    # r1 = a.run(r4, ICP=["eps", "n", "tau"], ISW=-1, DS=0.01)
    # r1 = r1 + a.run(r4, ICP=["eps", "n", "tau"], ISW=-1, DS=-0.01)  
    # rtot = rtot + r1

rtot = a.merge(rtot)

# a.plot(rtot, bifurcation_x="n", bifurcation_y="eps", bifurcation_z="tau", use_labels=False)

# --------------------
# r1 = a.run(prob, ICP=["eps", "n", "tau"])
# r1 = r1 + a.run(prob, ICP=["eps", "n", "tau"], DS=-0.01)
# rtot = r1

# r2 = a.run(r1('LP1'), ICP=["tau", "eps", "n"], DS=-0.01)
# rtot = rtot + r2
# r1 = a.run(r2, ICP=["eps", "n", "tau"], DS=0.01)
# r1 = r1 + a.run(r2, ICP=["eps", "n", "tau"], DS=-0.01)
# rtot = rtot + r1

# r2 = a.run(r1('LP1'), ICP=["tau", "eps", "n"], DS=-0.01)

# rtot = a.merge(rtot)
# a.plot(rtot, bifurcation_x="n", bifurcation_y="eps", bifurcation_z="tau", use_labels=False)
# --------------------


# rtau = a.run(r1, ICP=["tau", "eps", "n"], NMX=500, RL0 = 2.1, DSMAX=0.2, DS=-0.1)
# r2 = a.run(prob, ICP=["eps", "n", "tau"], RL0 = 0, DSMAX=0.5, DS=0.1)

 #a.plot(rtot, bifurcation_x="N", bifurcation_y="eps", bifurcation_z="T", use_labels=False)
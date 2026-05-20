import auto as a

sysname  = "wrinkle"
unames   = dict( enumerate( ["W", "Wp", "Wpp", "Wppp", "x"], start=1) )
parnames = dict( enumerate( ["eps", "n", "tau"], start=1) )

prob = a.load(
	e = sysname,
    unames = unames,
    parnames = parnames,
    dat='wrinkle.dat',
    NDIM=   len(unames), IPS = 4, IRS = 0, ILP = 1,
    ICP =  ["eps", "n", "tau"],
    NTST=  100, NCOL= 4, IAD = 3, ISP = 3, ISW = 1, IPLT= 0, NBC= 6, NINT= 0,
    NMX=  100, RL0 = 0, RL1 = 1.e+8, A0 = -1.e+8, A1 = 1.e+8, 
    NPR=  10, MXBF = 10, IID = 2, ITMX= 8, ITNW= 5, NWTN= 3, JAC= 0,
    EPSL= 1.e-8, EPSU = 1.e-8, EPSS = 1.e-6,
    DS  =  0.01, DSMIN=  0.001, DSMAX= 0.1, IADS= 1,
    NPAR = 5, THL = {}, THU = {},
    )

# We first run the code backwards to detect the first fold point
r1 = a.run(prob, SP=['LP1'])

# Beginning from fold point, solve downwards in tau
rtau = a.run(r1, IRS='LP1', SP=['LP'], ICP=["tau", "eps", "n"], NMX=200, RL0=10.0, DS=-0.1, DSMAX=0.1, LAB=1)

# Now for each tau run, find the fold point to the right
nmat = [rtau["n"]]
epsmat = [rtau["eps"]]
taumat = [rtau["tau"]]
for j in range(1, len(rtau())):
	# print "Solution label", rtau["LAB"]
	solution = rtau(j)
	rtmp = a.run(a.load(solution, SP=['LP1'], RL0 = 0.0, DS=0.001, DSMAX=0.1, ICP=["eps", "n", "tau"]))
	if j == 1:
		nmat = [rtmp["n"]]
		epsmat = [rtmp["eps"]]
		taumat = [rtmp["tau"]] 
	else:
		nmat.append(rtmp["n"])
		epsmat.append(rtmp["eps"])
		taumat.append(rtmp["tau"])

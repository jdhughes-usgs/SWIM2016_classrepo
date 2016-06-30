import numpy as np
import flopy
import os
import config

def henryfun(model_ws, dmcoef, nlay, ncol):
    # model grid information
    #nlay = 20
    nrow = 1
    #ncol = 42
    top = 1.
    bot = 0.
    dx = 2. / ncol
    dz = (top - bot) / nlay
    botm = np.array([[z] * ncol for z in np.arange(1-dz, 1-(nlay + 1) * dz, -dz)]).reshape((nlay, nrow, ncol))
    delr = np.array([dx] * ncol)
    delr[-1] = 0.01
    delc = 1

    # temporal discretization
    perlen = 0.15
    nstp = 1
    tsmult = 1.0
    ttsmult = 1.1

    # Input variables
    hk = 864.
    vk = 864.
    porosity = 0.35
    ss = 0.0001
    #dmcoef = 1.629
    Qinflow = 5.702
    Csalt = 35.0001
    Cfresh = 0.
    densesalt = 1025.
    densefresh = 1000.
    denseslp = (densesalt - densefresh) / (Csalt - Cfresh)

    # Create the data needed to make the constant head package
    chddata = []
    for k in range(nlay):
        chddata.append([k, 0, ncol - 1, 1.0, 1.0])

    # Create the data needed to make the well package
    weldata = []
    wellQ = Qinflow / nlay
    for k in range(nlay):
        weldata.append([k, 0, 0, wellQ])

    # create the ssm data
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    ssmdata = []
    for rec in chddata:
        k = rec[0]
        i = rec[1]
        j = rec[2]
        ssmdata.append([k, i, j, Csalt, itype['CHD']])
    for rec in weldata:
        k = rec[0]
        i = rec[1]
        j = rec[2]
        ssmdata.append([k, i, j, Cfresh, itype['WEL']])

    ssmdata = {0: ssmdata}
    chddata = {0: chddata}
    weldata = {0: weldata}

    # create the flopy objects
    modelname = 'henry'
    m = flopy.seawat.Seawat(modelname, model_ws=model_ws, exe_name=config.swexe)

    # modflow packages
    dis = flopy.modflow.ModflowDis(m, nlay, nrow, ncol, nper=1, delr=delr,
                                   delc=delc, laycbd=0, top=top,
                                   botm=botm, perlen=perlen, nstp=nstp)
    bas = flopy.modflow.ModflowBas(m)
    lpf = flopy.modflow.ModflowLpf(m, hk=hk, vka=hk)
    chd = flopy.modflow.ModflowChd(m, stress_period_data=chddata)
    wel = flopy.modflow.ModflowWel(m, stress_period_data=weldata)
    pcg = flopy.modflow.ModflowPcg(m, hclose=1.e-8)
    oc = flopy.modflow.ModflowOc(m, stress_period_data={(0, 0): ['save head', 'save budget']},
                                 compact=True)

    # mt3d packages
    btn = flopy.mt3d.Mt3dBtn(m, nprs=-1, prsity=porosity, sconc=Csalt, ifmtcn=0,
                             chkmas=False, nprobs=10, nprmas=10, dt0=1.e-4, ttsmult=ttsmult)
    adv = flopy.mt3d.Mt3dAdv(m, mixelm=0)
    dsp = flopy.mt3d.Mt3dDsp(m, al=0., trpt=1., trpv=1., dmcoef=dmcoef)
    gcg = flopy.mt3d.Mt3dGcg(m, iter1=500, mxiter=1, isolve=1, cclose=1e-7)
    ssm = flopy.mt3d.Mt3dSsm(m, stress_period_data=ssmdata)

    # seawat packages
    vdf = flopy.seawat.SeawatVdf(m, iwtable=0, densemin=0, densemax=0,
                                 denseref=densefresh, denseslp=denseslp, firstdt=1e-3)

    # write the SEAWAT input files
    m.write_input()

    try:
        fname = os.path.join(model_ws, 'MT3D001.UCN')
        os.remove(fname)
    except:
        pass
    if os.path.isfile(fname):
        raise Exception('Cannot proceed UCN file still exists.')
    success, dummy = m.run_model(silent=True)
    if not success:
        raise Exception('Model did not run successfully')

    fname = os.path.join(model_ws, 'henry.cbc')
    budobj = flopy.utils.CellBudgetFile(fname)
    chdflows = budobj.get_data(text='CONSTANT HEAD')[0]
    chdflows.q

    #caluculate the sum of the negative flows
    chdinflow = 0.
    for value in chdflows.q:
        if value > 0.:
            chdinflow = chdinflow + value
    return chdinflow
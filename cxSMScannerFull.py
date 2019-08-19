from __future__ import print_function
from __future__ import division
import cxSMTransitionsPhysBasis as cxSM
import random as rd
import argparse
import time
from math import sqrt


parser = argparse.ArgumentParser(prog='cxSMScanner',description='Scan the parameters in cxSM for EWPT.')
parser.add_argument('-n',dest='npoints',default=-1,type=int)
parser.add_argument('-i',dest='inputfile',default='cxSMPoints.dat')
parser.add_argument('-o',dest='outfile',default='cxSMScanned_Try01.dat')
args = parser.parse_args()

NPOINTS = args.npoints
NID = 0
NTRIED = 0
model = cxSM.cxSM(vS=100.0,mHH=500.0,mHA=500.0,theta=0.1,a1=-100.0,Min=[])



# Please Check whether following Vc and Tc functions are correct, changing them if necessary
def GetVLc(trans):
    if trans['crit_trans'] is None:
        return -1,-1
    else:
        return trans['crit_trans']['low_vev'][0],trans['crit_trans']['low_vev'][1]

def GetVHc(trans):
    if trans['crit_trans'] is None:
        return -1,-1
    else:
        return trans['crit_trans']['high_vev'][0],trans['crit_trans']['high_vev'][1]

def GetTc(trans):
    if trans['crit_trans'] is None:
        return -1
    else:
        return trans['crit_trans']['Tcrit']

def GetTypeC(trans):
    if trans['crit_trans'] is None:
        return -1
    else:
        return trans['crit_trans']['trantype']

def GetTypeN(trans):
    return trans['trantype']

def GetVLn(trans):
    return trans['low_vev'][0],trans['low_vev'][1]

def GetVHn(trans):
    return trans['high_vev'][0],trans['high_vev'][1]

def GetTn(trans):
    return trans['Tnuc']

# timeout = time.time() + 60*10

with open(args.outfile,'w') as outfile:
    with open(args.inputfile,'r') as inputfile:
        for line in inputfile:
            part = line.split()
            ID = int(part[0])
            vSin = float(part[1])
            mHHin = float(part[2])
            mHAin = float(part[3])
            thetain = float(part[4])
            a1in = float(part[5])
            # vSin = rd.uniform(50.0,2000.0)
            # mHHin = rd.uniform(200.0,5000.0)
            # mHAin = rd.uniform(200.0,5000.0)
            # thetain = rd.uniform(0.01,0.30)
            # a1in = -rd.uniform(10,500)**3

            "%f  %f  %f  %f  %f  "%(vSin,mHHin,mHAin,thetain,a1in)
            NTRIED+=1
            model.SetParameters(VS=vSin,mHH=mHHin,mHA=mHAin,theta=thetain,a1=a1in)
            print(model.PrintPhysicalParameters() + model.PrintPotentialParameters())

            if model.Stability() and model.Unitarity():
                try:
                    phases=model.getPhases()
                # except NotImplementedError as err:
                    # print 'Error in remove redundant phases'
                except:
                    print('Error in getPhases()')
                    continue
                if len(phases) == 0:
                    continue
                print(phases)
                TcTrans=model.calcTcTrans()
                TransString = 'TC  '
                print('=========================Tc len: ',len(TcTrans),'================================')
                # if len(Trans) == 0:
                #     continue
                TransString += '%d  '%(len(TcTrans))
                for i in range(len(TcTrans)):
                    TypeC=TcTrans[i]['trantype']
                    vhHC,vsHC=TcTrans[i]['high_vev'][0],TcTrans[i]['high_vev'][1]
                    vhLC,vsLC=TcTrans[i]['low_vev'][0],TcTrans[i]['low_vev'][1]
                    TC=TcTrans[i]['Tcrit']
                    TransString += '%d  %f  %f  %f  %f  %f  '%(TypeC,vhLC,vsLC,vhHC,vsHC,TC)
                try:
                    Trans = model.findAllTransitions()
                except:
                    print('Some Error, Skip this point')
                    continue
                TransString += 'TN  '
                print('=========================Tn len: ',len(Trans),'================================')
                if len(Trans) == 0:
                    continue
                TransString += '%d  '%(len(Trans))
                for i in range(len(Trans)):
                    TypeN=GetTypeN(Trans[i])
                    vhLN,vsLN=GetVLn(Trans[i])
                    vhHN,vsHN=GetVHn(Trans[i])
                    TN=GetTn(Trans[i])
                    TypeC=GetTypeC(Trans[i])
                    vhLC,vsLC=GetVLc(Trans[i])
                    vhHC,vsHC=GetVHc(Trans[i])
                    TC=GetTc(Trans[i])
                    alpha=Trans[i]['alpha_GW']
                    betaHn=Trans[i]['betaHn_GW']
                    TransString += '%d  %.10f  %.10f  %.10f  %.10f  %.10f  %d  %.10f  %.10f  %.10f  %.10f  %.10f  %.10f  %.10f  '%(TypeN,vhLN,vsLN,vhHN,vsHN,TN,TypeC,vhLC,vsLC,vhHC,vsHC,TC,alpha,betaHn)
                # if len(Trans) == 1:
                #     TransString = '1  %d  %f  %f  %d  %f  %f  '%(Trans[0]['trantype'],Trans[0]['action'],Trans[0]['Tnuc'],0,0,0)
                # elif len(Trans) == 2:
                #     TransString = '2  %d  %f  %f  %d  %f  %f  '%(Trans[0]['trantype'],Trans[0]['action'],Trans[0]['Tnuc'],Trans[1]['trantype'],Trans[1]['action'],Trans[1]['Tnuc'])

                # test = 0
                # if test == 5 or time.time() > timeout:
                #     break
                # test = test - 1
                    
                outfile.write("%d  "%ID+model.PrintPhysicalParameters()+model.PrintPotentialParameters()+TransString+'\n')
                outfile.flush()
                NID += 1
                if NPOINTS > 0 and NID >= NPOINTS:
                    print("Tried ",NTRIED," points, got ",NID," points.")
                    break

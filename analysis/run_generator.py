#!/usr/bin/env python3

import safari_input   # Makes the input files
import os             # Makes directories, etc
import numpy          # Used for frange
import argparse       # parses command line arguments

def frange(start, end, step):
    return numpy.arange(start, end, step)

def generate(basename, tstart, tend, pstart, pend, estart, eend, safio, res=0.025,\
             estep=1, pstep=1, tstep=1, TSA=None):
    
    if TSA is not None:
        TSA = float(TSA)
    
    rundir = os.path.join('.','runs')
    try:
        os.makedirs(rundir)
    except:
        pass
    rundir = os.path.join(rundir,basename)
    try:
        os.makedirs(rundir)
    except:
        pass
        
    for theta in frange(tstart, tend + tstep, tstep):
        for phi in frange(pstart, pend + pstep, pstep):
            for e0 in frange(estart, eend + estep, estep):
                edir = os.path.join(rundir,str(e0))
                try:
                    os.makedirs(edir)
                except:
                    pass
                safio.E0 = e0
                safio.THETA0 = theta
                safio.PHI0 = phi
                
                if TSA is not None:
                    safio.DTECTPAR[0] = 180 - (TSA + safio.THETA0)
                
                # Set energy resolution to 1%
                safio.ESIZE = e0/100
                filename = str(e0)+'_'+str(theta)\
                                  +'_'+str(phi)+'.input'
                filename = os.path.join(edir,filename)
                safio.save(file=filename)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file names")

    parser.add_argument("--ei", help="Inital Energy")
    parser.add_argument("--ef", help="Final Energy")
    parser.add_argument("--de", help="Energy Step")

    parser.add_argument("--ti", help="Inital Theta")
    parser.add_argument("--tf", help="Final Theta")
    parser.add_argument("--dt", help="Theta Step")

    parser.add_argument("--pi", help="Inital Phi")
    parser.add_argument("--pf", help="Final Phi")
    parser.add_argument("--dp", help="Phi Step")

    parser.add_argument("--tsa", help="Total Scattering Angle")

    parser.add_argument("-e", help="Energy (acts for both ei and ef)")
    parser.add_argument("-t", help="Theta (acts for both ti and tf)")
    parser.add_argument("-p", help="Phi (acts for both pi and pf)")

    parser.add_argument("-f","--face", help="Surface Face")

    parser.add_argument("-i","--input", help="Input File to use as a Template")

    parser.add_argument("-n","--number", help="Trajectory number (Monte Carlo) or Thermal Iterations (Adaptive Grid)")
    
    parser.add_argument("--montecarlo", help="Sets to Monte Carlo mode", action='store_true')
    parser.add_argument("--adaptivegrid", help="Sets to Adaptive Grid mode, argument of the value to use for depth!")
    
    args = parser.parse_args()

    basename = args.output

    estart = args.ei
    eend = args.ef
    estep = args.de

    tstart = args.ti
    tend = args.tf
    tstep = args.dt

    pstart = args.pi
    pend = args.pf
    pstep = args.dp

    face_line = args.face

    tsa = args.tsa
    
    template = 'template.input'
    
    if args.input is not None:
        template = args.input
        if not '.input' in template:
            template = template + '.input'

    if args.e is not None:
        estart = args.e
        eend = args.e

    if args.t is not None:
        tstart = args.t
        tend = args.t

    if args.p is not None:
        pstart = args.p
        pend = args.p

    if basename is None:
        basename = input("Output Base Name? ")
    
    if estart is None:
        estart = input("Energy Start? ")
    if eend is None:
        eend = input("Energy End? ")

    if estep is None:
        estep = '1'
    if args.de is None and estart!=eend:
        estep = input("Energy Step? ")
    
    if tstart is None:
        tstart = input("Theta Start? ")
    if tend is None:
        tend = input("Theta End? ")
    if tstep is None:
        tstep = '1'
    if args.dt is None and tstart!=tend:
        tstep = input("Theta Step? ")
    
    if pstart is None:
        pstart = input("Phi Start? ")
    if pend is None:
        pend = input("Phi End? ")
    if pstep is None:
        pstep = '1'
    if args.dp is None and pstart!=pend:
        pstep = input("Phi Step? ")

    if face_line is None:
        face_line = input("Surface Face (h k l)? ")
    
    face = None
    if face_line != '':
        face = safari_input.parseLine(face_line)
        
    safio = safari_input.SafariInput(template)
    
    if args.number is not None:
        safio.NUMCHA = int(args.number)
    
    if args.adaptivegrid is not None:
        safio.setAdaptiveGrid(int(adaptivegrid))
    elif args.montecarlo:
        safio.setMonteCarlo(True)
    
    if face is not None:
        safio.face = face
        safio.load_crystal = False
    else:
        safio.load_crystal = True
    
    generate(basename, float(tstart), float(tend), float(pstart),\
                       float(pend), float(estart), float(eend), safio,\
                       estep=float(estep), pstep=float(pstep), tstep=float(tstep), TSA=tsa)
